### Clear environment
rm(list = ls())

### Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(cluster)
library(igraph)

### Load clustered object
seurat_obj <- readRDS("pbmc_clustered.rds")

### Subset cluster 0
cluster0_cells <- WhichCells(seurat_obj, idents = "0")
cluster0_cells <- cluster0_cells[seurat_obj@meta.data[cluster0_cells, "cell_label"] != "cytotoxic_t"]
cluster0_cells <- cluster0_cells[seurat_obj@meta.data[cluster0_cells, "cell_label"] != "cd14_monocytes"]
expr <- GetAssayData(seurat_obj[, cluster0_cells], slot = "counts")

### Sweep parameters
nfeatures_list <- c(seq(40, 500, by = 20), seq(600, 4000, by = 100))
resolutions <- seq(0.02, 1.0, by = 0.02)
fixed_ndims <- 7  # Use fixed number of PCs

### Store results
all_results <- list()

for (nfeatures in nfeatures_list) {
  message(sprintf("Processing nfeatures: %d | ndims: %d", nfeatures, fixed_ndims))
  
  object <- CreateSeuratObject(expr)
  object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE,
                        variable.features.n = nfeatures)
  object <- RunPCA(object, npcs = max(fixed_ndims, 30), verbose = FALSE)
  object <- RunUMAP(object, reduction = "pca", dims = 1:fixed_ndims, verbose = FALSE)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:fixed_ndims, verbose = FALSE)
  
  for (res in resolutions) {
    object <- FindClusters(object, resolution = res, verbose = FALSE)
    clusters <- as.integer(Idents(object))
    
    if (length(unique(clusters)) > 1) {
      # Silhouette score
      pca_mat <- Embeddings(object, "pca")
      dmat <- as.matrix(dist(pca_mat[, 1:fixed_ndims]))
      sil <- silhouette(clusters, dmat)
      sil_score <- mean(sil[, 3])
      
      # Modularity score from SNN graph
      snn_graph <- object@graphs$SCT_snn
      ig <- graph_from_adjacency_matrix(as.matrix(snn_graph), mode = "undirected", weighted = TRUE)
      modularity_score <- modularity(ig, clusters)
    } else {
      sil_score <- NA
      modularity_score <- NA
    }
    
    all_results[[length(all_results) + 1]] <- data.frame(
      group = "cluster0",
      nfeatures = nfeatures,
      ndims = fixed_ndims,
      resolution = res,
      silhouette = sil_score,
      modularity = modularity_score
    )
  }
}

### Combine results and save
silhouette_df <- do.call(rbind, all_results)
write.csv(silhouette_df, "cluster0_silhouette_scores_fixed_ndim7.csv", row.names = FALSE)
