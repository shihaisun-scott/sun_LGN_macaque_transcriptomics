### Clear environment
rm(list = ls())

### Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(cluster)
library(mgcv)
library(clue)
library(Matrix)

### Load clustered Seurat object
obj <- readRDS("pbmc_clustered.rds")

### Subset cluster 0
cluster0_cells <- WhichCells(obj, idents = "0")
cluster0_cells <- cluster0_cells[obj@meta.data[cluster0_cells, "cell_label"] != "cytotoxic_t"]
cluster0_cells <- cluster0_cells[obj@meta.data[cluster0_cells, "cell_label"] != "cd14_monocytes"]

expr <- GetAssayData(obj[, cluster0_cells], slot = "counts")
cell_labels_cluster0 <- obj$cell_label[cluster0_cells]

# Variance explained for increasing feature steps
feature_steps <- seq(500, 4000, by = 500)
pc_variances <- list()

for (f in feature_steps) {
  cat("Processing", f, "features\n")
  
  obj_tmp <- CreateSeuratObject(expr)
  obj_tmp$cell_label <- cell_labels_cluster0
  obj_tmp <- SCTransform(obj_tmp, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = f)
  obj_tmp <- RunPCA(obj_tmp, npcs = 30, verbose = FALSE)
  
  sdev <- obj_tmp[["pca"]]@stdev
  var_explained <- sdev^2 / sum(sdev^2)
  
  pc_variances[[as.character(f)]] <- var_explained
}

# Construct dataframe from results
pc_df <- do.call(rbind, lapply(names(pc_variances), function(f) {
  data.frame(
    PC = 1:length(pc_variances[[f]]),
    Variance = pc_variances[[f]],
    Features = as.numeric(f)
  )
}))

# Plot separate elbow plots
ggplot(pc_df, aes(x = PC, y = Variance)) +
  geom_line() +
  facet_wrap(~ Features, scales = "free_y") +
  labs(title = "Elbow Plots for Different Numbers of Features",
       x = "Principal Component",
       y = "Proportion of Variance Explained")

# Combined plot with overlapping lines
pc_df$Features <- factor(pc_df$Features)

overlap_plot <- ggplot(pc_df, aes(x = PC, y = Variance, group = Features, color = Features)) +
  geom_line(alpha = 0.7) +
  geom_vline(xintercept = 7, linetype = "dashed", color = "red") +
  labs(title = "Overlapping Elbow Plots for Varying Feature Counts",
       x = "Principal Component",
       y = "Proportion of Variance Explained") +
  theme_minimal()

# Save plot
pdf("1_pc_detection_overlap.pdf", height = 5, width = 10)
print(overlap_plot)
dev.off()

# Display plot in RStudio
overlap_plot
