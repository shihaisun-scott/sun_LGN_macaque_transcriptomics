# perform clustering and save seurat objects
# set directory to this script's folder

rm(list = ls())
gc()
set.seed(10)


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
# resolutions <- c(M = 0.2, P = 0.26, K = 0.12) # the original optimized res
resolutions <- c(M = 0.3, P = 0.3, K = 0.26) # produces the same clusters but higher res
ndims <- c(M = 4, P = 4, K = 4)


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
  
  # --- Sensitivity Analysis ---
  match_stats <- data.frame()
  original_n_clusters <- length(unique(Idents(object)))
  skip_match_count <- 0
  
  
  for (iter in 1:1000) {
    
    # define the random sampling here
    p_vals <- c(p_cells = runif(1, 0.7, 0.9), p_nfeatures = runif(1, 0.8, 1.2), p_ndims = sample(c(-1, 0, 1), 1))
    
    # 1. Random downsampling percentage between 70–90%
    p <- p_vals["p_cells"]
    sampled_cells <- sample(colnames(object), size = floor(p * ncol(object)))
    
    # 2. Downsample object
    sampled_obj <- subset(object, cells = sampled_cells)
    
    # 3. Randomize nfeatures (80–120%)
    alt_nfeatures <- round(nfeatures * p_vals["p_nfeatures"])

    # 4. Randomize ndims ±1
    alt_ndim <- ndim + p_vals["p_ndims"]

    # Reprocess
    sampled_obj[["RNA"]] <- split(sampled_obj[["RNA"]], f = sampled_obj$donor)
    sampled_obj <- SCTransform(sampled_obj, verbose = FALSE, return.only.var.genes = FALSE,
                               variable.features.n = alt_nfeatures)
    sampled_obj <- RunPCA(sampled_obj, npcs = 30, verbose = FALSE)
    sampled_obj <- IntegrateLayers(sampled_obj, method = HarmonyIntegration,
                                   orig.reduction = "pca", new.reduction = "harmony",
                                   verbose = FALSE, assay = "SCT")
    sampled_obj[["RNA"]] <- JoinLayers(sampled_obj[["RNA"]])
    sampled_obj <- RunUMAP(sampled_obj, reduction = "harmony", dims = 1:alt_ndim, verbose = FALSE)
    sampled_obj <- FindNeighbors(sampled_obj, reduction = "harmony", dims = 1:alt_ndim, verbose = FALSE)
    sampled_obj <- FindClusters(sampled_obj, resolution = res, verbose = FALSE)
    
    # Match clusters only if cluster number matches original
    sampled_n_clusters <- length(unique(Idents(sampled_obj)))
    
    if (sampled_n_clusters == original_n_clusters) {
      clusters <- Idents(sampled_obj)
      original_subset <- object$cluster_id[names(clusters)]
      contingency <- table(clusters, original_subset)
      contingency_mat <- as.matrix(contingency)
      padded <- matrix(0, nrow = max(dim(contingency_mat)), ncol = max(dim(contingency_mat)))
      padded[1:nrow(contingency_mat), 1:ncol(contingency_mat)] <- contingency_mat
      assignment <- solve_LSAP(max(padded) - padded)
      matched <- sum(padded[cbind(1:nrow(padded), assignment)])
      match_percent <- matched / length(clusters) * 100
      
      # cat(sprintf("Iter %02d | %.0f%% match | %.0f%% cells | %d feats | %d dims\n",
      #             iter, match_percent, p * 100, alt_nfeatures, alt_ndim))
    } else {
      match_percent <- NA
      skip_match_count <- skip_match_count + 1
      # cat(sprintf("Iter %02d | Skipped (cluster mismatch: %d vs %d)\n",
      #             iter, sampled_n_clusters, original_n_clusters))
    }
    
    match_stats <- rbind(match_stats, data.frame(iteration = iter,
                                                 percent_cells = round(p * 100),
                                                 match = match_percent,
                                                 features_used = alt_nfeatures,
                                                 ndims_used = alt_ndim))
  }
  
  # Report stats
  summary_stats <- match_stats %>%
    summarise(mean_match = mean(match, na.rm = TRUE), sd_match = sd(match, na.rm = TRUE))
  
  print(summary_stats)
  cat(sprintf("Total skipped due to cluster mismatch: %d\n", skip_match_count))
  
  write.csv(match_stats, sprintf("analysis_output/3b_sensitivity_analysis_%s.csv", nam), row.names = FALSE)
}
  