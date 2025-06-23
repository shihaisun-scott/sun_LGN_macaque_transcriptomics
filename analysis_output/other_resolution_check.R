rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(mgcv)
library(tidyr)

load("data/data.rda")
load("data/precluster_info.rda")


subdat=datlist[,precluster$sample_name]
plot_objects=list()
subcells=list()
subcells[["M"]] <- subset(metalist$sample_name, precluster$cluster_label == "M")
subcells[["P"]] <- subset(metalist$sample_name, precluster$cluster_label == "P")
subcells[["K"]] <- subset(metalist$sample_name, grepl("^K", precluster$cluster_label))

df <- read.csv("analysis_output/silhouette_scores_sweep_2.csv", stringsAsFactors = FALSE) 

# custom settings
feature <- c(M = 2800, P = 3200, K = 1100)
resolutions <- c(M = 0.2, P = 0.26, K = 0.12)
ndims <- c(M = 4, P = 4, K = 4)

for (nam in c("M","P","K")){
  fe <- feature[nam]
  res <- resolutions[nam]
  ndim <- ndims[nam]
  sil <- df[df$group == nam & df$nfeatures == fe & df$resolution == res & df$ndims == ndim, "silhouette"]
  
  same_res <- na.omit(df[df$group == nam & df$nfeatures == fe & df$ndims == ndim & df$silhouette == sil, "resolution"])
  
  
  # also check if percentage of cells are the same with the changed resolution
  keepcells <- subcells[[nam]]
  object <- CreateSeuratObject(counts = subdat[, keepcells], meta.data = metalist[keepcells, ])
  object[["RNA"]] <- split(object[["RNA"]], f = object$donor)
  object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = fe)
  object <- RunPCA(object, npcs = 30, verbose = FALSE)
  object <- IntegrateLayers(object, method = HarmonyIntegration, orig.reduction = "pca", 
                            new.reduction = "harmony", verbose = FALSE, assay = "SCT")
  object[["RNA"]] <- JoinLayers(object[["RNA"]])
  object <- RunUMAP(object, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
  object <- FindNeighbors(object, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
  object <- FindClusters(object, resolution = res, verbose = FALSE)
  
  # second set with different res
  sampled_obj <- object
  sampled_obj <- FindClusters(sampled_obj, resolution = max(same_res), verbose = FALSE)
  
  object$cluster_id <- Idents(object)
  sampled_obj$cluster_id <- Idents(sampled_obj)
  
  sampled_n_clusters <- length(unique(Idents(sampled_obj)))
  
  clusters <- Idents(sampled_obj)
  original_subset <- object$cluster_id[names(clusters)]
  contingency <- table(clusters, original_subset)
  contingency_mat <- as.matrix(contingency)
  padded <- matrix(0, nrow = max(dim(contingency_mat)), ncol = max(dim(contingency_mat)))
  padded[1:nrow(contingency_mat), 1:ncol(contingency_mat)] <- contingency_mat
  assignment <- solve_LSAP(max(padded) - padded)
  matched <- sum(padded[cbind(1:nrow(padded), assignment)])
  match_percent <- matched / length(clusters) * 100
  
  
  cat(sprintf("%s | max res = %f | match = %.0f%% \n", nam, max(same_res), match_percent))
  
  
  
}