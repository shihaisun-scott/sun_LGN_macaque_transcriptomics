rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(cluster)
library(mgcv)
library(tidyr)

# -------------------------------
# Load data
# -------------------------------
load("D:/Partners HealthCare Dropbox/Shi Sun/Research/Pezaris/Documents/Manuscripts/transcriptomics/Rscripts/data/data.rda")
load("D:/Partners HealthCare Dropbox/Shi Sun/Research/Pezaris/Documents/Manuscripts/transcriptomics/Rscripts/data/precluster_info.rda")
load("D:/Partners HealthCare Dropbox/Shi Sun/Matlab scripts/Pezaris/transcriptomy/macaque_lgn_2021_txn/Rscripts/data/sct_newclusters.rda")

# Parameters
n_dims <- c(4, 5, 6, 7, 8)
n_k_val <- c(3, 5, 7, 10)

# Cell subsets
subdat <- datlist[, precluster$sample_name]
subcells <- list(
  M = subset(metalist$sample_name, precluster$cluster_label == "M"),
  P = subset(metalist$sample_name, precluster$cluster_label == "P"),
  K = subset(metalist$sample_name, grepl("^K", precluster$cluster_label))
)

# Methods
# methods_to_plot <- c("max_silhouette", "smoothed_gam")
methods_to_plot <- "smoothed_gam"

# Loop over groups and methods to collect plots for each
for (group in names(subcells)) {
  for (method in methods_to_plot) {
    
    all_plots <- list()
    
    for (ndim in n_dims) {
      for (k_val in n_k_val) {
        
        df <- read.csv("analysis_output/silhouette_scores_sweep_2.csv", stringsAsFactors = FALSE) %>%
          filter(!is.na(silhouette), resolution > 0.1, ndims == ndim, nfeatures > 1000, nfeatures < 3900) %>%
          filter(nfeatures %% 10 != 5)
        
        if (nrow(df) == 0) next
        
        # Select best param
        best_max <- df %>%
          group_by(group) %>%
          slice_max(silhouette, n = 1, with_ties = FALSE) %>%
          mutate(method = "max_silhouette")
        
        best_gam <- df %>%
          group_by(group) %>%
          group_modify(~{
            model <- gam(silhouette ~ s(nfeatures, k = k_val) + s(resolution, k = k_val), data = .x)
            .x$predicted <- predict(model)
            .x %>% slice_max(predicted, n = 1)
          }) %>%
          mutate(method = "smoothed_gam")
        
        all_methods_clean <- bind_rows(best_max, best_gam)
        
        # Get row for this group & method
        param_row <- all_methods_clean %>%
          filter(group == !!group, method == !!method)
        
        if (nrow(param_row) == 0) next
        
        keepcells <- subcells[[group]]
        nfeatures <- param_row$nfeatures[[1]]
        ndim_param <- param_row$ndims[[1]]
        res <- param_row$resolution[[1]]
        
        # Seurat pipeline
        object <- CreateSeuratObject(counts = subdat[, keepcells], meta.data = newclusters[keepcells, ])
        object[["RNA"]] <- split(object[["RNA"]], f = object$donor)
        object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE,
                              variable.features.n = nfeatures)
        object <- RunPCA(object, npcs = 30, verbose = FALSE)
        object <- IntegrateLayers(object, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony",
                                  verbose = FALSE, assay = "SCT")
        object[["RNA"]] <- JoinLayers(object[["RNA"]])
        object <- RunUMAP(object, reduction = "harmony", dims = 1:ndim_param, verbose = FALSE)
        object <- FindNeighbors(object, reduction = "harmony", dims = 1:ndim_param, verbose = FALSE)
        object <- FindClusters(object, resolution = res, verbose = FALSE)
        object$cluster_id <- Idents(object)
        
        # Silhouette
        sil_obj <- silhouette(as.integer(object$cluster_id),
                              dist(Embeddings(object, "harmony")[, 1:ndim_param]))
        sil_score <- round(mean(sil_obj[, 3]), 3)
        
        # Sensitivity
        n_boot <- 10
        boot_scores <- numeric(n_boot)
        original_clusters <- Idents(object)
        
        for (b in 1:n_boot) {
          set.seed(100 + b)
          downsampled <- sample(colnames(object), size = floor(0.8 * ncol(object)))
          obj_boot <- subset(object, cells = downsampled)
          obj_boot <- FindNeighbors(obj_boot, reduction = "harmony", dims = 1:ndim_param, verbose = FALSE)
          obj_boot <- FindClusters(obj_boot, resolution = res, verbose = FALSE)
          overlap <- intersect(downsampled, names(original_clusters))
          match <- sum(as.integer(original_clusters[overlap]) == as.integer(Idents(obj_boot)[overlap]))
          boot_scores[b] <- match / length(overlap)
        }
        
        sensitivity_score <- round(mean(boot_scores), 2)
        
        param_label <- paste0("Dims: ", ndim_param,
                              ", Feat: ", nfeatures,
                              ", Res: ", round(res, 2),
                              "\nSil: ", sil_score,
                              ", Sens: ", sensitivity_score)
        
        # Cluster UMAP
        plot_cluster <- DimPlot(object, group.by = "cluster_id", reduction = "umap", alpha = 0.5) +
          labs(x = param_label, y = NULL) +
          ggtitle(paste0("Cluster - ndim ", ndim, ", k ", k_val)) +
          theme_void() +
          theme(aspect.ratio = 1,
                plot.title = element_text(hjust = 0.5, size = 9),
                axis.title.x = element_text(size = 7),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                legend.position = "none")
        
        all_plots[[paste(ndim, k_val, sep = "_")]] <- plot_cluster
      }
    }
    
    # Save 1 PDF per group-method
    if (length(all_plots) > 0) {
      final_grid <- plot_grid(plotlist = all_plots, ncol = 5)
      savename <- paste0("UMAP_", group, "_", method, "_all_ndim_kval.pdf")
      pdf(savename, width = 13, height = 8)
      print(final_grid)
      dev.off()
      cat("Saved:", savename, "\n")
    }
  }
}

