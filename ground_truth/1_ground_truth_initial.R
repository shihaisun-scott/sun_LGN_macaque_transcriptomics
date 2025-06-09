# plot umap with standard seurat pipeline

rm(list = ls())


library(Seurat)
library(cowplot)

expr<-readRDS("pbmc_900.rds")
cell_label<-read.table("cell_labels.txt",header = T)
cell_label<-factor(cell_label$cell_type_labels,levels = c("regulatory_t","naive_t","cd34","memory_t","cd56_nk","cytotoxic_t","cd14_monocytes","b_cells","naive_cytotoxic"))

obj<-CreateSeuratObject(expr)
obj$cell_label<-cell_label
obj<-NormalizeData(obj,verbose=FALSE)
obj<-FindVariableFeatures(obj,verbose=FALSE)
obj<-ScaleData(obj,verbose=FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj<-RunUMAP(obj,dims = 1:30,verbose=FALSE)

umap1 <- DimPlot(obj,reduction = "umap",group.by = "cell_label",shuffle = T)

# clustering and save
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)

umap2 <- DimPlot(obj,reduction = "umap",group.by = "ident",shuffle = T)


# save both umaps as figure
combined_plot <- plot_grid(umap1, umap2, labels = "AUTO", ncol = 2)
ggsave(filename = paste0("standard_umap.pdf"), plot = combined_plot, width = 16, height = 6)


# Save the final Seurat object for analysis down the line
saveRDS(obj, file = "pbmc_clustered.rds")
