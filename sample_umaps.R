# perform clustering and save seurat objects
# set directory to this script's folder

rm(list = ls())


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)


load("D:/Partners HealthCare Dropbox/Shi Sun/Research/Pezaris/Documents/Manuscripts/transcriptomics/Rscripts/data/data.rda")
load("D:/Partners HealthCare Dropbox/Shi Sun/Research/Pezaris/Documents/Manuscripts/transcriptomics/Rscripts/data/precluster_info.rda")

# load the old clusters
load("D:/Partners HealthCare Dropbox/Shi Sun/Matlab scripts/Pezaris/transcriptomy/macaque_lgn_2021_txn/Rscripts/data/sct_newclusters.rda")


subdat=datlist[,precluster$sample_name]
# subdat = datlist[setdiff(rownames(datlist), "NSUN6"), precluster$sample_name]


plot_objects=list()
subcells=list()
subcells[["M"]] <- subset(metalist$sample_name, precluster$cluster_label == "M")
subcells[["P"]] <- subset(metalist$sample_name, precluster$cluster_label == "P")
subcells[["K"]] <- subset(metalist$sample_name, grepl("^K", precluster$cluster_label))


# custom settings

# GAM surface; also pretty good results
# feature <- c(M = 3600, P = 700, K = 600)
# ndims <- c(M = 10, P = 10, K = 10)
# resolutions <-  c(M = 0.34, P = 0.2, K = 0.2)

# this is good - 05/20
feature <- c(M = 2400, P = 2300, K = 1100)
ndims <- c(M = 10, P = 10, K = 10)
resolutions <-  c(M = 0.4, P = 0.28, K = 0.38)


plot_list <- list()
plot_list2 <- list()


seurat_objects=list()


# loop between each set of cells
for (nam in names(subcells)) {
  count = 0
  keepcells=subcells[[nam]]
  nfeatures = feature[[nam]]
  ndim = ndims[[nam]]
  res = resolutions[[nam]]
  
  object <- CreateSeuratObject(counts = subdat[,keepcells], meta.data = newclusters[keepcells,])
  object[["RNA"]] <- split(object[["RNA"]], f = object$donor)
  object <- SCTransform(object = object, verbose = FALSE, return.only.var.genes = FALSE,
                        variable.features.n = nfeatures)
  object <- RunPCA(object = object, npcs = 30, verbose = FALSE)
  object <- IntegrateLayers(object = object, method = HarmonyIntegration,
                            orig.reduction = "pca", new.reduction = "harmony",
                            verbose = FALSE, assay = "SCT")
  #object[["RNA"]] <- JoinLayers(object[["RNA"]])
  #object[["SCT"]] <- JoinLayers(object[["SCT"]])
  
  object <- RunUMAP(object, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
  object <- FindNeighbors(object, reduction = "harmony", dims = 1:ndim, verbose = FALSE)
  object <- FindClusters(object, resolution = res, verbose = FALSE)
  object$cluster_id <- Idents(object)
  
  
  # plot UMAP
  dim_plot <- DimPlot(object, group.by = "ident", alpha = 0.5) + 
    ggtitle(nam) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme_void() 
  
  plot_list[[nam]] <- dim_plot
  
  dim_plot <- DimPlot(object, group.by = "cluster_label", alpha = 0.5) + 
    ggtitle(nam) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme_void() 
  
  plot_list2[[nam]] <- dim_plot
  
}


library(cowplot)

# Adjust layout (2 rows × 2 columns, change as needed)
n_cols <- 3
n_rows <- ceiling(length(plot_list) / n_cols)

pdfname <- "2_sample_umaps_max2.pdf"
pdf(pdfname, height = 3 * n_rows, width = 3 * n_cols)

# Combine all plots into one grid
combined_list <- c(plot_list, plot_list2)
combined_plot <- plot_grid(plotlist = combined_list, ncol = n_cols)

# Save grid to PDF
print(combined_plot)

dev.off()


combined_plot



