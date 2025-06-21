# global search for optimum nfeatures and clustering resolution based on silhouette scores

rm(list = ls())

# Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(cluster)
library(igraph)

# load data from 1_data_collection
load("data/data.rda")
load("data/precluster_info.rda")

# categorize groups
subcells <- list(
  M = subset(metalist$sample_name, precluster$cluster_label == "M"),
  P = subset(metalist$sample_name, precluster$cluster_label == "P"),
  K = subset(metalist$sample_name, grepl("^K", precluster$cluster_label))
)

subdat <- datlist[, precluster$sample_name]


# Parameter ranges
nfeatures_list <- seq(500, 5000, by = 100)
ndims_list <- c(4, 5, 6, 7, 8, 10)
resolutions <- seq(0.02, 1.0, by = 0.02)

# Store all silhouette results
all_results <- list()

# Run sweep
for (nam in names(subcells)) {
  keepcells <- subcells[[nam]]
  
  for (nfeatures in nfeatures_list) {
    for (ndims in ndims_list) {
      message(sprintf("Processing %s | nfeatures: %d | ndims: %d", nam, nfeatures, ndims))
      
      object <- CreateSeuratObject(counts = subdat[, keepcells], meta.data = metalist[keepcells, ])
      object[["RNA"]] <- split(object[["RNA"]], f = object$donor)
      object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE,
                            variable.features.n = nfeatures)
      object <- RunPCA(object, npcs = 30, verbose = FALSE)
      object <- IntegrateLayers(object, method = HarmonyIntegration,
                                orig.reduction = "pca", new.reduction = "harmony",
                                verbose = FALSE, assay = "SCT")
      object[["RNA"]] <- JoinLayers(object[["RNA"]])
      object <- RunUMAP(object, reduction = "harmony", dims = 1:ndims, verbose = FALSE)
      object <- FindNeighbors(object, reduction = "harmony", dims = 1:ndims, verbose = FALSE)
      
      for (res in resolutions) {
        object <- FindClusters(object, resolution = res, verbose = FALSE)
        object$cluster_res <- Idents(object)
        
        clusters <- as.integer(object$cluster_res)
        
        if (length(unique(clusters)) > 1) {
          
          # silhouette score
          pca_mat <- Embeddings(object, "harmony")[, 1:ndims]
          dist_mat <- dist(pca_mat)
          sil <- silhouette(clusters, dist_mat)
          sil_score <- mean(sil[, 3])
          
          # Modularity
          snn_graph <- object@graphs$SCT_snn
          ig <- graph_from_adjacency_matrix(as.matrix(snn_graph), mode = "undirected", weighted = TRUE)
          modularity_score <- modularity(ig, clusters)
        } else {
          sil_score <- NA
          modularity_score <- NA
        }
        
        all_results[[length(all_results) + 1]] <- data.frame(
          group = nam,
          nfeatures = nfeatures,
          ndims = ndims,
          resolution = res,
          silhouette = sil_score,
          modularity = modularity_score
          
        )
      }
    }
  }
}

# Combine results
silhouette_df <- do.call(rbind, all_results)

# Find best parameter set per group
best_params <- silhouette_df %>%
  group_by(group, nfeatures, ndims, resolution) %>%
  summarize(silhouette = mean(silhouette, na.rm = TRUE), .groups = "drop") %>%
  group_by(group) %>%
  slice_max(silhouette, n = 1, with_ties = FALSE)

print("Best parameter set per group:")
print(best_params)

write.csv(silhouette_df, "analysis_output/silhouette_scores_sweep_2.csv", row.names = FALSE)

