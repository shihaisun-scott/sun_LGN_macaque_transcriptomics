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

### Load clustered Seurat object
obj <- readRDS("pbmc_clustered.rds")

### Subset cluster 0
cluster0_cells <- WhichCells(obj, idents = "0")
cluster0_cells <- cluster0_cells[obj@meta.data[cluster0_cells, "cell_label"] != "cytotoxic_t"]
cluster0_cells <- cluster0_cells[obj@meta.data[cluster0_cells, "cell_label"] != "cd14_monocytes"]

expr <- GetAssayData(obj[, cluster0_cells], slot = "counts")
cell_labels_cluster0 <- obj$cell_label[cluster0_cells]

### Set parameters
nfeatures = 1500
ndim = 7
res = 0.24

nfeatures = 3900
ndim = 5
res = 0.44

# standard : can get 85.71% with these parameters
# nfeatures = 2000
# ndim = 20
# res = 0.8

### Create Seurat object with selected features
object <- CreateSeuratObject(expr)
object$cell_label <- cell_labels_cluster0
object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = nfeatures)
object <- RunPCA(object, npcs = 20, verbose = FALSE)

# elbow plot to check PCs
ElbowPlot(object)

### UMAP + Clustering
object <- RunUMAP(object, dims = 1:ndim, verbose = FALSE)
object <- FindNeighbors(object, dims = 1:ndim, verbose = FALSE)
object <- FindClusters(object, resolution = res, verbose = FALSE)



### also calculate the percentage matching
# clusters <- Idents(object)
# labels <- object$cell_label
# confusion <- table(clusters, labels)
# correct_matches <- sum(apply(confusion, 1, max))
# total_cells <- length(clusters)
# percent_match <- correct_matches / total_cells * 100
# cat(sprintf(">> Percentage match between clusters and true labels: %.2f%%\n", percent_match))

# Create contingency table
clusters <- Idents(object)
labels <- object$cell_label
confusion <- table(clusters, labels)
conf_mat <- as.matrix(confusion)
n <- max(nrow(conf_mat), ncol(conf_mat))
square_mat <- matrix(0, n, n)
square_mat[1:nrow(conf_mat), 1:ncol(conf_mat)] <- conf_mat
cost_mat <- max(square_mat) - square_mat
assignment <- solve_LSAP(cost_mat)
correct_matches <- sum(square_mat[cbind(1:n, assignment)])
total_cells <- length(clusters)
percent_match <- correct_matches / total_cells * 100
cat(sprintf(">> Optimally matched percentage: %.2f%%\n", percent_match))



# plot umaps and save with percentage match
umap1 <- DimPlot(object, group.by = "ident", label = TRUE, repel = TRUE) + 
  ggtitle(sprintf("%g features, %g PCs, %g res: Match = %.2f%%", nfeatures, ndim, res, percent_match))

umap2 <- DimPlot(object, group.by = "cell_label", label = TRUE, repel = TRUE) +
  ggtitle("Ground Truth Cell Types")

combined_plot <- plot_grid(umap1, umap2, labels = "AUTO", ncol = 2)
ggsave(filename = paste0("matching_umap.pdf"), plot = combined_plot, width = 16, height = 6)


