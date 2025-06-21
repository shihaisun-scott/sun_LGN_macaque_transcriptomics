# perform clustering and save seurat objects
# set directory to this script's folder

rm(list = ls())


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(clue) 


load("data/data.rda")
load("data/precluster_info.rda")


subdat=datlist[,precluster$sample_name]
plot_objects=list()
subcells=list()
subcells[["M"]] <- subset(metalist$sample_name, precluster$cluster_label == "M")
subcells[["P"]] <- subset(metalist$sample_name, precluster$cluster_label == "P")
subcells[["K"]] <- subset(metalist$sample_name, grepl("^K", precluster$cluster_label))


# custom settings
feature <- c(M = 2800, P = 3200, K = 1100)
resolutions <- c(M = 0.2, P = 0.26, K = 0.24)
ndims <- c(M = 4, P = 4, K = 7)


newclusters = precluster
newclusters$cluster_id <- c(matrix(0, 1, length(precluster$sample_name)))
cluster_count = 0;

plot_list <- list()
seurat_objects=list()

subdat=datlist[][,precluster$sample_name]

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
  
  # sensitivity analysis
  #downsample_percents <- c(0.7, 0.8, 0.9)
  downsample_percents <- c(0.8)
  match_stats <- data.frame()
  
  for (p in downsample_percents) {
    for (iter in 1:100) {

      set.seed(123 + iter)
      
      # downsample cells
      sampled_cells <- sample(colnames(object), size = floor(p * ncol(object)))
      
      # stratified
      #sampled_cells <- unlist(lapply(split(colnames(object), object$cluster_id), function(cells_in_cluster) {
      #  sample(cells_in_cluster, size = floor(p * length(cells_in_cluster)))
      #}))
      
      sampled_obj <- subset(object, cells = sampled_cells)
      
      # Vary number of features for sensitivity
      alt_nfeatures <- sample(c(nfeatures, round(nfeatures * 1.2), round(nfeatures * 0.8)), 1)
      #alt_nfeatures <- nfeatures
      
      sampled_obj[["RNA"]] <- split(sampled_obj[["RNA"]], f = sampled_obj$donor)
      sampled_obj <- SCTransform(sampled_obj, verbose = FALSE, return.only.var.genes = FALSE,
                                 variable.features.n = alt_nfeatures)
      sampled_obj <- RunPCA(sampled_obj, npcs = 30, verbose = FALSE)
      sampled_obj <- IntegrateLayers(sampled_obj, method = HarmonyIntegration,
                                     orig.reduction = "pca", new.reduction = "harmony",
                                     verbose = FALSE, assay = "SCT")
      sampled_obj[["RNA"]] <- JoinLayers(sampled_obj[["RNA"]])
      sampled_obj <- RunUMAP(sampled_obj, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
      sampled_obj <- FindNeighbors(sampled_obj, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
      sampled_obj <- FindClusters(sampled_obj, resolution = res, verbose = FALSE)
      
      # Compare clusters
      clusters <- Idents(sampled_obj)
      original_subset <- object$cluster_id[names(clusters)]
      contingency <- table(clusters, original_subset)
      contingency_mat <- as.matrix(contingency)
      padded <- matrix(0, nrow = max(dim(contingency_mat)), ncol = max(dim(contingency_mat)))
      padded[1:nrow(contingency_mat), 1:ncol(contingency_mat)] <- contingency_mat
      assignment <- solve_LSAP(max(padded) - padded)
      matched <- sum(padded[cbind(1:nrow(padded), assignment)])
      match_percent <- matched / length(clusters) * 100
      
      match_stats <- rbind(match_stats, data.frame(iteration = iter, percent = p * 100,
                                                   match = match_percent, features_used = alt_nfeatures))
      
      cat(sprintf(" %.0f%% \n", match_percent))
      
    }
  }
  
  # Report stats
  summary_stats <- match_stats %>%
    group_by(percent) %>%
    summarise(mean_match = mean(match), sd_match = sd(match))
  print(summary_stats)
  
  write.csv(match_stats, sprintf("analysis_output/sensitivity_analysis_%s.csv", nam), row.names = FALSE)
  
}



