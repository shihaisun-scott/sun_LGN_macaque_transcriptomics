# figures for parameter optimization - elbow, gam heatmap, umap

rm(list = ls())


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(Matrix)
library(mgcv)
library(tidyr)

load("data/data.rda")
load("data/precluster_info.rda")

dat <- datlist[, precluster$sample_name]

groups <- list(
  M = subset(metalist$sample_name, precluster$cluster_label == "M"),
  P = subset(metalist$sample_name, precluster$cluster_label == "P"),
  K = subset(metalist$sample_name, grepl("^K", precluster$cluster_label))
)

umap_params <- list(
  M = list(nfeatures = 2800, ndims = 4, resolution = 0.2),
  P = list(nfeatures = 3200, ndims = 4, resolution = 0.26),
  K = list(nfeatures = 1100, ndims = 7, resolution = 0.24)
)

ndims_map <- c(M = 4, P = 4, K = 7)

sil_df <- read.csv("analysis_output/silhouette_scores_sweep_2.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(silhouette), resolution > 0.1, nfeatures > 1000, nfeatures < 3900)

# Store plot objects
elbow_data_list <- list()
heatmap_list <- list()
umap_list <- list()

# loop over groups
for (grp in names(groups)) {
  message("Processing: ", grp)
  cells <- groups[[grp]]
  param <- umap_params[[grp]]
  
  # elbow
  feature_steps <- seq(1000, 4000, by = 500)
  pc_variances <- list()
  
  for (f in feature_steps) {
    expr <- dat[, cells]
    meta <- metalist[cells, , drop = FALSE]
    obj <- CreateSeuratObject(expr, meta.data = meta)
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$donor)
    obj <- SCTransform(obj, variable.features.n = f, verbose = FALSE, return.only.var.genes = FALSE)
    obj <- RunPCA(obj, npcs = 15, verbose = FALSE)
    sdev <- obj[["pca"]]@stdev
    pc_variances[[as.character(f)]] <- sdev^2 / sum(sdev^2)
  }
  
  pc_df <- do.call(rbind, lapply(names(pc_variances), function(f) {
    data.frame(PC = 1:length(pc_variances[[f]]),
               Variance = pc_variances[[f]],
               Features = as.numeric(f),
               Group = grp)
  }))
  elbow_data_list[[grp]] <- pc_df
  
  # gam heatmap
  sil_df_group <- sil_df %>%
    filter(group == grp, ndims == ndims_map[[grp]]) %>%
    mutate(resolution_round = round(resolution, 2))
  
  gam_model <- gam(silhouette ~ s(nfeatures, k = 10) + s(resolution_round, k = 10), data = sil_df_group)
  
  full_grid <- expand.grid(nfeatures = unique(sil_df_group$nfeatures),
                           resolution_round = unique(sil_df_group$resolution_round))
  full_grid$silhouette_pred <- predict(gam_model, newdata = full_grid)
  
  # Mark the max predicted silhouette point
  best_point <- full_grid %>% slice_max(silhouette_pred, n = 1)
  
  heatmap <- ggplot(full_grid, aes(x = factor(nfeatures), y = factor(resolution_round))) +
    geom_tile(aes(fill = silhouette_pred)) +
    geom_point(data = best_point, aes(x = factor(nfeatures), y = factor(resolution_round)),
               color = "red", shape = 21, size = 3, stroke = 1) +
    geom_text(data = best_point,
              aes(label = round(silhouette_pred, 3)),
              size = 2.5, vjust = -1.2, color = "black") +
    scale_fill_viridis_c(option = "D", limits = c(0, NA)) +
    labs(title = paste("GAM Grid search -", grp),
         x = "Features", y = "Resolution") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  heatmap_list[[grp]] <- heatmap
  
  # umap
  expr <- dat[, cells]
  meta <- metalist[cells, , drop = FALSE]
  obj <- CreateSeuratObject(expr, meta.data = meta)
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$donor)
  obj <- SCTransform(obj, variable.features.n = param$nfeatures, verbose = FALSE, return.only.var.genes = FALSE)
  obj <- RunPCA(obj, npcs = 15, verbose = FALSE)
  obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", assay = "SCT", verbose = FALSE)
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:param$ndims)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:param$ndims)
  obj <- FindClusters(obj, resolution = param$resolution)
  obj$cluster_id <- Idents(obj)
  
  # silhouette score
  sil_row <- sil_df_group %>%
    filter(nfeatures == param$nfeatures,
           resolution_round == round(param$resolution, 2)) %>%
    slice_max(silhouette, n = 1)
  sil_score <- ifelse(nrow(sil_row) > 0, round(sil_row$silhouette, 3), NA)
  
  umap_title <- paste0(grp, " - UMAP\n",
                       "Features: ", param$nfeatures,
                       ", Dims: ", param$ndims,
                       ", Res: ", param$resolution,
                       ", Silhouette: ", sil_score)
  
  umap_plot <- DimPlot(obj, group.by = "cluster_id", reduction = "umap") +
    ggtitle(umap_title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 9))
  
  umap_list[[grp]] <- umap_plot
}

# combine elbow plots for overlapping lines
all_pc_df <- bind_rows(elbow_data_list)
chosen_ndims <- data.frame(Group = c("M", "P", "K"),
                           ndims = c(5, 7, 7))

elbow_plot <- ggplot(all_pc_df, aes(x = PC, y = Variance)) +
  geom_line() +
  geom_vline(data = chosen_ndims, aes(xintercept = ndims), color = "red", linetype = "dashed") +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  labs(title = "Elbow Plots with Chosen Dimensions",
       x = "Principal Component", y = "Variance Explained") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

elbow_list <- list()
chosen_ndims <- c(M = 4, P = 4, K = 7)

for (grp in names(elbow_data_list)) {
  df <- elbow_data_list[[grp]]
  df$Features <- factor(df$Features)
  
  elbow_list[[grp]] <- ggplot(df, aes(x = PC, y = Variance, color = Features, group = Features)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    geom_vline(xintercept = chosen_ndims[[grp]], linetype = "dashed", color = "red") +
    labs(title = paste("Elbow Plot -", grp),
         x = "Principal Component", y = "Variance Explained",
         color = "N Features") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7))
}

# now combine all figs
final_plot_grid <- plot_grid(
  elbow_list$M, heatmap_list$M, umap_list$M,
  elbow_list$P, heatmap_list$P, umap_list$P,
  elbow_list$K, heatmap_list$K, umap_list$K,
  ncol = 3,
  label_size = 10
)

# save
pdf("analysis_output/2d_post_search.pdf", width = 16, height = 12)
print(final_plot_grid)
dev.off()
