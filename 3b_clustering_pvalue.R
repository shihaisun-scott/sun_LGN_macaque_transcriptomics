# shuffle transcripts and calculate silhouette scores for p-value

rm(list = ls())


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(cluster)



load("data/data.rda")
load("data/precluster_info.rda")


subdat=datlist[,precluster$sample_name]

plot_objects=list()
subcells=list()
subcells[["M"]] <- subset(metalist$sample_name, precluster$cluster_label == "M")
subcells[["P"]] <- subset(metalist$sample_name, precluster$cluster_label == "P")
subcells[["K"]] <- subset(metalist$sample_name, grepl("^K", precluster$cluster_label))


# custom settings
feature <- c(M = 2600, P = 2900, K = 1100)
resolutions <- c(M = 0.22, P = 0.24, K = 0.26)
ndims <- c(M = 4, P = 5, K = 7)

newclusters = precluster
newclusters$cluster_id <- c(matrix(0, 1, length(precluster$sample_name)))
cluster_count = 0;

plot_list <- list()
seurat_objects=list()

# loop between each set of cells
for (nam in names(subcells)) {
  keepcells <- subcells[[nam]]
  nfeatures <- feature[[nam]]
  ndim <- ndims[[nam]]
  res <- resolutions[[nam]]
  
  object <- CreateSeuratObject(counts = subdat[, keepcells], meta.data = metalist[keepcells, ])
  object[["RNA"]] <- split(object[["RNA"]], f = object$donor)
  object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = nfeatures)
  object <- RunPCA(object, npcs = 30, verbose = FALSE)
  object <- IntegrateLayers(object, method = HarmonyIntegration, orig.reduction = "pca", 
                            new.reduction = "harmony", verbose = FALSE, assay = "SCT")
  object[["RNA"]] <- JoinLayers(object[["RNA"]])
  object <- RunUMAP(object, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
  object <- FindNeighbors(object, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
  object <- FindClusters(object, resolution = res, verbose = FALSE)
  object$cluster_id <- Idents(object)
  
  # Save plot
  plot_list[[nam]] <- DimPlot(object, group.by = "ident", alpha = 0.5) +
    ggtitle(nam) + 
    theme_void() +
    theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  
  
  
  # Calculate distance matrix and silhouette score
  real_dist <- dist(Embeddings(object, reduction = "harmony")[, 1:ndim])
  real_clusters <- as.integer(object$cluster_id)
  sil <- silhouette(real_clusters, real_dist)
  observed_silhouette_score <- mean(sil[, 3])  # average silhouette width

  
  # Cache raw data once
  raw_counts <- as.matrix(subdat[, keepcells])
  n_geneperm_iter <- 100
  perm_silhouette_scores <- numeric(n_geneperm_iter)
  
  for (i in 1:n_geneperm_iter) {
    set.seed(1+i)
    # Shuffle expression of each gene across cells
    permuted_counts <- apply(raw_counts, 1, sample)
    permuted_counts <- matrix(permuted_counts, nrow = nrow(raw_counts), byrow = FALSE)
    rownames(permuted_counts) <- rownames(raw_counts)
    colnames(permuted_counts) <- colnames(raw_counts)
    
    # Create Seurat object from permuted data
    perm_obj <- CreateSeuratObject(counts = permuted_counts, meta.data = metalist[colnames(permuted_counts), ])
    perm_obj[["RNA"]] <- split(perm_obj[["RNA"]], f = object$donor)
    perm_obj <- SCTransform(perm_obj, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = nfeatures)
    perm_obj <- RunPCA(perm_obj, npcs = 30, verbose = FALSE)
    perm_obj <- IntegrateLayers(perm_obj, method = HarmonyIntegration, orig.reduction = "pca", 
                                new.reduction = "harmony", verbose = FALSE, assay = "SCT")
    perm_obj[["RNA"]] <- JoinLayers(perm_obj[["RNA"]])
    perm_obj <- RunUMAP(perm_obj, dims = 1:ndim, verbose = FALSE)
    perm_obj <- FindNeighbors(perm_obj, dims = 1:ndim, verbose = FALSE)
    perm_obj <- FindClusters(perm_obj, resolution = res, verbose = FALSE)
    
    perm_clusters <- as.integer(Idents(perm_obj))
    
    if (length(unique(perm_clusters)) > 1) {
      perm_dist <- dist(Embeddings(perm_obj, reduction = "harmony")[, 1:ndim])
      sil <- silhouette(perm_clusters, perm_dist)
      perm_silhouette_scores[i] <- mean(sil[, 3])
    } else {
      perm_silhouette_scores[i] <- NA
    }
  }
  # Remove NA values
  perm_silhouette_scores <- na.omit(perm_silhouette_scores)
  
  # Compute p-value
  p_val_silhouette <- mean(perm_silhouette_scores >= observed_silhouette_score)
  
  cat(sprintf("Silhouette p-value for %s: %.4f\n", nam, p_val_silhouette))
  
}




