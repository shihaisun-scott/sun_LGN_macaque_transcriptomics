# perform clustering and save seurat objects
# set directory to this script's folder

rm(list = ls())


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)


load("data/data.rda")
load("data/precluster_info.rda")


subdat=datlist[,precluster$sample_name]
# subdat = datlist[rownames(datlist) != "NSUN6", precluster$sample_name] # do not include nsun6

plot_objects=list()
subcells=list()
subcells[["M"]] <- subset(metalist$sample_name, precluster$cluster_label == "M")
subcells[["P"]] <- subset(metalist$sample_name, precluster$cluster_label == "P")
subcells[["K"]] <- subset(metalist$sample_name, grepl("^K", precluster$cluster_label))
subcells[["MPK"]] <- subset(metalist$sample_name, 
                            precluster$cluster_label == "M" | 
                              precluster$cluster_label == "P" | 
                              grepl("^K", precluster$cluster_label))


# custom settings
feature <- c(M = 2800, P = 3200, K = 1100, MPK = 3000)
resolutions <- c(M = 0.2, P = 0.26, K = 0.24, MPK = 0.2)
ndimensions <- c(M = 4, P = 4, K = 7, MPK = 15)


newclusters = precluster
newclusters$cluster_id <- c(matrix(0, 1, length(precluster$sample_name)))
cluster_count = 0;

plot_list <- list()
seurat_objects=list()


# loop between each set of cells
for (nam in names(subcells)) {
  count = 0
  keepcells=subcells[[nam]]
  nfeatures = feature[[nam]]
  ndims = ndimensions[[nam]]
  cluster_res = resolutions[[nam]]
  
  object <- CreateSeuratObject(counts = subdat[,keepcells], meta.data = metalist[keepcells,])
  object[["RNA"]] <- split(object[["RNA"]], f = object$donor)
  object <- SCTransform(object = object, verbose = FALSE, return.only.var.genes = FALSE,
                        variable.features.n = nfeatures)
  object <- RunPCA(object = object, npcs = 30, verbose = FALSE)
  object <- IntegrateLayers(object = object, method = HarmonyIntegration,
                            orig.reduction = "pca", new.reduction = "harmony",
                            verbose = FALSE, assay = "SCT")
  object[["RNA"]] <- JoinLayers(object[["RNA"]])
  object <- RunUMAP(object, reduction = "harmony", dims = 1:ndims)
  object <- FindNeighbors(object, reduction = "harmony", dims = 1:ndims)
  object <- FindClusters(object, resolution = cluster_res)
  object$cluster_id <- Idents(object)
  
  
  # plot UMAP
  dim_plot <- DimPlot(object, group.by = "ident", alpha = 0.5) + 
    ggtitle(nam) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme_void() 
  
  plot_list[[nam]] <- dim_plot
  
  
  # save new clustering
  if (nam %in% c("M", "P", "K")){
    unique(object@active.ident)
    cluster_results <- list()
    
    for (cluster_id in unique(object@active.ident)){
      cluster_name <- paste0(nam, as.numeric(cluster_id)+1) # +1 to avoid the 0
      newclusters[names(object@active.ident)[object@active.ident == cluster_id],"cluster_label"] <- cluster_name
      
      cluster_count = cluster_count + 1
      newclusters[names(object@active.ident)[object@active.ident == cluster_id],"cluster_id"] <- cluster_count
    }
  }
  seurat_objects[[nam]] = object
}

# save new clusters
save(newclusters,file="data/sct_newclusters.rda")

# save seurat objects
save(seurat_objects,file="data/sct_seurat_objects.rda")

# save umap plots
# pdfname <- paste0("analysis_output/sct_sample_umaps.pdf")
# pdf(pdfname, height = 12, width = 18)
# grid.arrange(grobs = plot_list)
# dev.off()

umaps <- plot_grid(plotlist = plot_list, ncol = 2)

pdfname <- paste0("analysis_output/3a_clustering_umaps.pdf")
pdf(pdfname, height = 6, width = 8)
print(umaps)
dev.off()





