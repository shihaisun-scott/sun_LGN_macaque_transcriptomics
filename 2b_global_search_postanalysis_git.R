# post analysis of global search -- determine the best parameters after GAM fitting

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(mgcv)
library(tidyr)

# load and filter
df <- read.csv("analysis_output/silhouette_scores_sweep_2.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(silhouette), resolution > 0.1, ndims == 4, nfeatures > 1000, nfeatures < 3900)

k_val <- 10 # df for GAM fitting
savename= "2c_global_search_postanalysis.pdf"

df <- df %>%
  filter(nfeatures %% 10 != 5)

# remove identical silhouette scores
#df <- df %>%
#  group_by(group, nfeatures, ndims, silhouette) %>%
#  slice_max(resolution, n = 1, with_ties = FALSE) %>%
#  ungroup()

# max silhouette score
best_max <- df %>%
  group_by(group) %>%
  slice_max(silhouette, n = 1, with_ties = FALSE) %>%
  mutate(method = "max_silhouette")


# smoothed gam surface
best_gam <- df %>%
  group_by(group) %>%
  group_modify(~{
    model <- gam(silhouette ~ s(nfeatures, k = k_val) + s(resolution, k = k_val), data = .x)
    .x$predicted <- predict(model)
    .x %>% slice_max(predicted, n = 1)
  }) %>%
  mutate(method = "smoothed_gam")


# combine
all_methods <- bind_rows(best_max, best_threshold, best_gam)

# clean some columns
all_methods_clean <- all_methods %>%
  select(where(~ !is.list(.)))

# save
print(all_methods_clean)
write.csv(all_methods_clean, "analysis_output/optimal_parameters_summary.csv", row.names = FALSE)


## Visualizae
# create heat map
df <- df %>%
  mutate(resolution_round = round(resolution, 2))

# fit GAM
df_filled <- df %>%
  group_by(group) %>%
  group_modify(~{
    data_group <- .x
    
    # Fit GAM model only on non-missing silhouette values
    model <- gam(silhouette ~ s(nfeatures, k = k_val) + s(resolution_round, k = k_val), data = data_group)
    
    # Create full grid of (nfeatures, resolution_round)
    full_grid <- expand.grid(
      nfeatures = unique(data_group$nfeatures),
      resolution_round = unique(data_group$resolution_round)
    )
    
    # Predict silhouette scores for all combinations
    full_grid$silhouette_pred <- predict(model, newdata = full_grid)
    
    # Merge predictions back with actual values
    full_grid <- full_grid %>%
      left_join(data_group, by = c("nfeatures", "resolution_round")) %>%
      mutate(group = unique(data_group$group)) %>%
      mutate(silhouette_combined = ifelse(is.na(silhouette), silhouette_pred, silhouette))
    
    return(full_grid)
  }) %>%
  ungroup()

# Ensure resolution is rounded the same way
all_methods_clean <- all_methods_clean %>%
  mutate(resolution_round = round(resolution, 2))

# Choose a label to display 
all_methods_clean <- all_methods_clean %>%
  mutate(label = method) 


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



### plot umaps here
load("data/data.rda")
load("data/precluster_info.rda")

# Load optimal parameter results
all_methods_clean <- read.csv("analysis_output/optimal_parameters_summary.csv", stringsAsFactors = FALSE)

subdat <- datlist[, precluster$sample_name]
subcells <- list()
subcells[["M"]] <- subset(metalist$sample_name, precluster$cluster_label == "M")
subcells[["P"]] <- subset(metalist$sample_name, precluster$cluster_label == "P")
subcells[["K"]] <- subset(metalist$sample_name, grepl("^K", precluster$cluster_label))

# methods to visualize
methods_to_plot <- c("max_silhouette", "smoothed_gam")

# Plot storage
plot_list <- list()

# umap loop
for (method in methods_to_plot) {
  for (nam in names(subcells)) {
    
    param_row <- all_methods_clean %>%
      filter(group == nam, method == !!method)
    
    if (nrow(param_row) == 0) {
      warning(paste("No parameters found for", nam, method))
      next
    }
    
    keepcells <- subcells[[nam]]
    nfeatures <- param_row$nfeatures[[1]]
    ndim <- param_row$ndims[[1]]
    res <- param_row$resolution[[1]]
    
    object <- CreateSeuratObject(counts = subdat[, keepcells], meta.data = newclusters[keepcells, ])
    object[["RNA"]] <- split(object[["RNA"]], f = object$donor)
    object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE,
                          variable.features.n = nfeatures)
    object <- RunPCA(object, npcs = 30, verbose = FALSE)
    object <- IntegrateLayers(object, method = HarmonyIntegration,
                              orig.reduction = "pca", new.reduction = "harmony",
                              verbose = FALSE, assay = "SCT")
    object[["RNA"]] <- JoinLayers(object[["RNA"]])
    object <- RunUMAP(object, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
    object <- FindNeighbors(object, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
    object <- FindClusters(object, resolution = res, verbose = FALSE)
    object$cluster_id <- Idents(object)
    
    param_label <- paste0("Dims: ", ndim, 
                          ", Features: ", nfeatures, 
                          ", Res: ", round(res, 2))
    
    # --- UMAP 1: Seurat cluster
    plot_cluster <- DimPlot(object, group.by = "cluster_id", alpha = 0.5) +
      ggtitle(paste(nam, "-", method)) +
      labs(x = param_label, y = NULL) +
      theme_void() +
      theme(aspect.ratio = 1,
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.title = element_text(hjust = 0.5, size = 10),
            legend.position = "none")
    
    # Store plots
    plot_list[[paste(nam, method, "cluster", sep = "_")]] <- plot_cluster
  }
}

umap_grid <- plot_grid(plotlist = plot_list, ncol = 3)

# Combine: heatmap on left, UMAPs on right
final_plot <- plot_grid(
  heatmap,  # left side
  umap_grid,
  ncol = 2,
  rel_widths = c(1.5, 2.5)  # adjust to preference
)

# Save as PDF
pdf(savename, width = 12, height = 8)
print(final_plot)
dev.off()


final_plot
