# -------------------------------
# Load libraries
# -------------------------------
library(dplyr)
library(ggplot2)
library(mgcv)
library(tidyr)
library(Seurat)
library(SeuratWrappers)
library(cowplot)

# -------------------------------
# Load silhouette results
# -------------------------------
df <- read.csv("cluster0_silhouette_scores_sweep.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(silhouette), resolution >= 0.2, ndims == 5, nfeatures > 499, nfeatures < 4000)

k_val = 5 # df for GAM fitting

savename <- "umap_methods_pc5.pdf"

# Add dummy group column if not present (for heatmap facet)
df$group <- "Cluster 0"

# -------------------------------
# Define methods
# -------------------------------
best_max <- df %>% slice_max(silhouette, n = 1, with_ties = FALSE) %>% mutate(method = "max_silhouette")

top_band <- df %>%
  mutate(threshold = quantile(silhouette, 0.95, na.rm = TRUE)) %>%
  filter(silhouette >= threshold) %>%
  arrange(desc(nfeatures), desc(ndims), resolution) %>%
  slice_head(n = 1) %>%
  mutate(method = "top_5_percent_band")

best_gam <- df %>%
  group_modify(~{
    model <- gam(silhouette ~ s(nfeatures, k = k_val) + s(resolution, k = k_val), data = .x)
    .x$predicted <- predict(model)
    .x %>% slice_max(predicted, n = 1)
  }) %>% mutate(method = "smoothed_gam")

best_threshold <- df %>%
  mutate(thresh = 0.95 * max(silhouette, na.rm = TRUE)) %>%
  filter(silhouette >= thresh) %>%
  arrange(desc(nfeatures), desc(ndims), resolution) %>%
  slice_head(n = 1) %>%
  mutate(method = "threshold_largest")

# Combine selected methods
all_methods <- bind_rows(best_max, best_threshold, best_gam)
all_methods_clean <- all_methods %>% select(where(~ !is.list(.)))

# Save to CSV
write.csv(all_methods_clean, "cluster0_parameters_summary.csv", row.names = FALSE)

# -------------------------------
# GAM Heatmap
# -------------------------------
df <- df %>% mutate(resolution_round = round(resolution, 2))
all_methods_clean <- all_methods_clean %>% mutate(resolution_round = round(resolution, 2), group = "Cluster 0")

df_filled <- df %>%
  group_by(group) %>%
  group_modify(~{
    model <- gam(silhouette ~ s(nfeatures, k = k_val) + s(resolution_round, k = k_val), data = .x)
    grid <- expand.grid(
      nfeatures = unique(.x$nfeatures),
      resolution_round = unique(.x$resolution_round)
    )
    grid$silhouette_pred <- predict(model, newdata = grid)
    grid <- grid %>%
      left_join(.x, by = c("nfeatures", "resolution_round")) %>%
      mutate(group = unique(.x$group)) %>%
      mutate(silhouette_combined = ifelse(is.na(silhouette), silhouette_pred, silhouette))
    return(grid)
  }) %>% ungroup()

heatmap <- ggplot(df_filled, aes(x = factor(nfeatures), y = factor(resolution_round))) +
  geom_tile(aes(fill = silhouette_combined)) +
  geom_point(data = all_methods_clean, aes(x = factor(nfeatures), y = factor(resolution_round)),
             color = "red", size = 3, shape = 21, fill = "white", stroke = 1) +
  geom_text(data = all_methods_clean, aes(x = factor(nfeatures), y = factor(resolution_round), label = method),
            color = "black", size = 2.5, vjust = -1) +
  scale_fill_viridis_c(option = "D") +
  facet_wrap(~ group) +
  labs(
    title = "Smoothed Silhouette Score Heatmap with Selected Methods",
    x = "Number of Features",
    y = "Resolution (rounded)",
    fill = "Silhouette"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )

# -------------------------------
# UMAP Generation
# -------------------------------
obj <- readRDS("pbmc_clustered.rds")
cluster0_cells <- WhichCells(obj, idents = "0")
cluster0_cells <- cluster0_cells[!(obj@meta.data[cluster0_cells, "cell_label"] %in% c("cytotoxic_t", "cd14_monocytes"))]

expr <- GetAssayData(obj[, cluster0_cells], slot = "counts")
cell_labels_cluster0 <- obj$cell_label[cluster0_cells]
methods_to_plot <- c("max_silhouette", "threshold_largest", "smoothed_gam")

plot_list <- list()
plot_list2 <- list()

for (method in methods_to_plot) {
  param_row <- all_methods_clean %>% filter(method == !!method)
  if (nrow(param_row) == 0) next
  
  nfeatures <- param_row$nfeatures[[1]]
  ndim <- param_row$ndims[[1]]
  res <- param_row$resolution[[1]]
  
  object <- CreateSeuratObject(expr)
  object$cell_label <- cell_labels_cluster0
  object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = nfeatures)
  object <- RunPCA(object, npcs = 20, verbose = FALSE)
  object <- RunUMAP(object, dims = 1:ndim, verbose = FALSE)
  object <- FindNeighbors(object, dims = 1:ndim, verbose = FALSE)
  object <- FindClusters(object, resolution = res, verbose = FALSE)
  object$cluster_id <- Idents(object)
  
  param_label <- paste0("Dims: ", ndim, ", Features: ", nfeatures, ", Res: ", round(res, 2))
  
  plot_cluster <- DimPlot(object, group.by = "cluster_id", alpha = 0.5) +
    ggtitle(paste("Seurat -", method)) +
    labs(x = param_label, y = NULL) +
    theme_void() +
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.position = "none")
  
  plot_label <- DimPlot(object, group.by = "cell_label", alpha = 0.5) +
    ggtitle(paste("True -", method)) +
    labs(x = param_label, y = NULL) +
    theme_void() +
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.position = "none")
  
  plot_list[[paste(method, "cluster")]] <- plot_cluster
  plot_list2[[paste(method, "label")]] <- plot_label
}

# Combine all UMAPs vertically
umap_grid <- plot_grid(plotlist = plot_list, ncol = 1)
umap_grid2 <- plot_grid(plotlist = plot_list2, ncol = 1)

# Combine heatmap and UMAPs side-by-side
final_plot <- plot_grid(
  heatmap,
  umap_grid,
  umap_grid2,
  ncol = 3,
  rel_widths = c(1.5, 2.5)
)

# Save combined figure
pdf(savename, width = 16, height = 8)
print(final_plot)
dev.off()

print(all_methods)
