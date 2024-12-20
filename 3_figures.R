# for each set of cells, plot umap, heat map, and dot plot

rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
require(Matrix)
library(SeuratWrappers)
library(gridExtra)
library(openxlsx)

load("data/data.rda")
load("data/precluster_info.rda")
load("data/sct_newclusters.rda")
load("data/sct_seurat_objects.rda")



## plot UMAP with bakken's clustering and color coding
object <- seurat_objects[["MPK"]]
mpk_colors <- c("#DB2D43","#0048BA","#177245", "#177245")
names(mpk_colors) <- c("M","P","K1", "K2")

num_m = 2 # number of m clusters
num_p = 3 # number of p clusters
num_k = 5 # number of k clusters

all_colors <- c(colorRampPalette(colors = c("darkred", "lightpink"))(num_m),
                colorRampPalette(colors = c("#2E2D88", "lightblue"))(num_p),
                colorRampPalette(colors = c("#004B49","#ACE1AF"))(4))
names(all_colors) <- c("M1","M2","P1","P2","P3","K1","K2","K3","K4")
barplot(rep(1, 9), col = all_colors, border = NA, names.arg = names(all_colors)) # check colors



# manually give cluster id and cluster names
object@meta.data$cluster_label <- metalist$cluster_label[match(rownames(object@meta.data), metalist$sample_name)]
object@meta.data$cluster_id <- metalist$cluster_id[match(rownames(object@meta.data), metalist$sample_name)]
object$predefined_clusters_label <- metalist$cluster_label[match(rownames(object@meta.data), metalist$sample_name)]
Idents(object) <- "cluster_id"

# plot UMAP
umap_plot <- DimPlot(object, group.by = "cluster_label", alpha = 0.5, cols = mpk_colors) + 
  ggtitle("") + 
  theme(aspect.ratio = 1,
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

pdf("analysis_output/umap_bakken.pdf", height = 5, width = 5)
print(umap_plot)
dev.off()

# exp_colors <- c("blue", "white", "red")
exp_colors <- c("purple", "black", "yellow")



## plot UMAP and heatmaps of the following: M cells, P cells, K cells, MP cells, MPK cells, all cells
subdat=datlist[,newclusters$sample_name]
subdat=sweep(subdat,2,Matrix::colSums(subdat),"/")*10^6
plot_objects=list()
subcells=list()
subcells[["M"]] <- subset(newclusters$sample_name, grepl("^M", newclusters$cluster_label))
subcells[["P"]] <- subset(newclusters$sample_name, grepl("^P", newclusters$cluster_label))
subcells[["K"]] <- subset(newclusters$sample_name, grepl("^K", newclusters$cluster_label))


for (nam in names(subcells)) {
  count = 0
  keepcells=subcells[[nam]]
  object <- seurat_objects[[nam]]
  plot_list = list()
  
  # manually give cluster id and cluster names
  object@meta.data$cluster_label <- newclusters$cluster_label[match(rownames(object@meta.data), newclusters$sample_name)]
  object@meta.data$cluster_id <- newclusters$cluster_id[match(rownames(object@meta.data), newclusters$sample_name)]
  object$predefined_clusters_label <- newclusters$cluster_label[match(rownames(object@meta.data), newclusters$sample_name)]
  
  Idents(object) <- "cluster_id"
  
  
  
  # plot UMAP
  umap_plot <- DimPlot(object, reduction = "umap", group.by = "cluster_label", alpha = 0.75, cols= all_colors, pt.size	= 3) + 
    ggtitle("") + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  plot_list[[1]] <- umap_plot
  
  # get differential gene expressions
  object <- PrepSCTFindMarkers(object,  assay = "SCT", verbose = FALSE)
  object.markers <- FindAllMarkers(object, assay = "SCT",
                                   test.use = "wilcox",  # Use Wilcoxon rank sum test (you can choose other tests)
                                   min.pct = 0.1,  # Minimum percent of cells expressing the gene
                                   logfc.threshold = 0.25, only.pos = TRUE)  # Log fold change threshold
  
  # only keep significant markers
  significant_markers <- object.markers %>%
    filter(p_val_adj < 0.05 )%>%
    arrange(cluster) 
  
  top_genes_per_cluster <- significant_markers %>%
    slice_max(order_by = avg_log2FC, n = Inf) %>%
    group_by(cluster)
  
  top_10_genes_per_cluster <- significant_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = -p_val_adj, n = 10)
  
  existing_genes <- rownames(GetAssayData(object, assay = "SCT", slot = "scale.data"))
  
  top_10_genes_filtered <- significant_markers %>%
    filter(gene %in% existing_genes) %>%
    group_by(cluster) %>%
    slice_max(order_by = -p_val_adj, n = 10)
  
  # order by avg_log2fc for top 30
  # keep top 10 for heat map visualization
  top30 <- significant_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10)
  
  
  
  # plot heat map
  heat_map <- DoHeatmap(object,
                        features = top_10_genes_filtered$gene,
                        group.by = "cluster_label",  # Ensure this matches your metadata column name
                        label = TRUE,
                        group.colors = all_colors, 
                        draw.lines = TRUE,
                        angle = 0,
                        raster = FALSE)
  
  heat_map3 <- heat_map +
    scale_fill_gradientn(colors = c("purple", "black", "yellow"))
  
  heat_map2 <- heat_map +
    scale_fill_gradientn(colors = c("blue", "white", "red"))
  
  heat_map1 <- heat_map +
    scale_fill_gradientn(colors = c("white", "blue"))
  
  plot_list[[2]] <- heat_map2
  
  
  # also plot dot plot with specific genes list
  if (nam %in% "M"){
    specific_genes <- c("CDC73","ZDHHC21","SOX5","PTPRD", "PCDH9",
                        "STX1A","NOTCH2","ABCA7","SPTBN4","PTPRS")
  }else if (nam %in% "P"){
    specific_genes <- c("NDUFA4","COX7C","ATP5ME","CYCS","GAPDH","NEFL","TUBA1B","ACTG1","CFL1","DYNLL1","HSP90AA1",
                        "NOTCH2","NRXN1","NRG3","PCDH9",
                        "LDLRAD4","KDM4B","KAZN","ADARB2")
  }else if (nam %in% "K"){
    specific_genes <- c("IL1RAPL2","DAPK1","MEGF10","SEL1L3","ADAMTSL1",
                        "MYO16","GJB6","SNTG2","ANK1","TTN","CCDC69","TSN1","EGFR","PLD1",
                        "KCNT2","KCNQ2","KCNK2","KCNIP1","CACNA1E","CACNG4","GRID2","GRIN3A","HRH1","OPRK1","TRPM3","MGLL",
                        "TNS3","PRKCA","VWA5B1")
  }
  
  dot_plot <- DotPlot(object,features=specific_genes,
                      group.by = "cluster_label",
                      scale.min = 0,
                      scale.max = 100, 
                      assay = "SCT") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  dot_plot2 <- dot_plot +
    # scale_color_gradientn(colors = c("blue", "lightgrey", "red"), limits = c(-2, 2))
    scale_color_gradientn(colors = c("blue", "white", "red"))
  
  
  dot_plot1 <- dot_plot +
    scale_color_gradientn(colors = c("white", "blue"))
  
  dot_plot3 <- dot_plot +
    scale_color_gradientn(colors = c("purple", "white", "yellow"))

  # dot_plot <- dot_plot + coord_flip()
  
  plot_list[[3]] <- dot_plot2
  
  # violin plots
  VlnPlot(object,
    features = specific_genes,
    ncol = length(specific_genes),
    same.y.lims = TRUE,
    group.by = "cluster_label",
    stack = TRUE)
  
  
  
  # Run UMAP based on the PCA of the selected features
  sig_genes <- significant_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = -p_val_adj, n = 140)
  significant_genes <- unique(significant_markers$gene)
  seurat_object <- RunUMAP(object, features = significant_genes)
  plot_list[[4]] <- DimPlot(seurat_object, reduction = "umap", group.by = "cluster_label", alpha = 0.75, cols= all_colors, pt.size = 3) + 
    ggtitle("") + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  
  # run umap again with a fresh seurat object
  new_object <- subset(object, features = significant_genes)
  new_object[["RNA"]] <- split(new_object[["RNA"]], f = object$donor)
  new_object <- SCTransform(object = new_object, verbose = FALSE, return.only.var.genes = TRUE,
                            variable.features.n = length(significant_genes))
  new_object <- RunPCA(object = new_object, npcs = 30, verbose = FALSE)
  new_object <- IntegrateLayers(object = new_object, method = HarmonyIntegration,
                                orig.reduction = "pca", new.reduction = "harmony",
                                verbose = FALSE, assay = "SCT")
  new_object[["RNA"]] <- JoinLayers(new_object[["RNA"]])
  new_object <- RunUMAP(new_object, reduction = "harmony", dims = 1:5)
  umap2 <- DimPlot(new_object, reduction = "umap", group.by = "cluster_label", alpha = 0.75, cols= all_colors, pt.size = 3) + 
    ggtitle(paste0("ngenes = ", length(significant_genes))) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  plot_list[[5]] <- umap2
  
  
  # save cluster, number of cells, markers into summary table
  cell_counts_per_cluster <- object@meta.data %>%
    group_by(predefined_clusters_label) %>%
    summarise(num_cells = n())
  
  # Ensure matching formats: convert cluster column to match predefined_clusters_label format
  top_genes_per_cluster <- top_genes_per_cluster %>%
    mutate(cluster = object@meta.data$cluster_label[match(cluster, Idents(object))])
  
  # Merge the top genes with the cell counts
  summary_table <- top_genes_per_cluster %>%
    group_by(cluster) %>%
    summarize(top_genes = paste(gene, collapse = ", ")) %>%
    left_join(cell_counts_per_cluster, by = c("cluster" = "predefined_clusters_label"))
  table_grob <- tableGrob(summary_table)
  
  tablename <- paste0("analysis_output/", nam, ".xlsx")
  
  
  
  # save the p-values for each significant gene in excel
  wb <- createWorkbook()
  addWorksheet(wb, "Summary")
  writeData(wb, sheet = "Summary", summary_table) # add aummary table to the first page
  clusters <- unique(top_genes_per_cluster$cluster)
  for (cluster in clusters) {
    cluster_genes <- top_genes_per_cluster %>% filter(cluster == !!cluster)
    addWorksheet(wb, paste0("Cluster_", cluster))
    writeData(wb, sheet = paste0("Cluster_", cluster), cluster_genes)
  }
  saveWorkbook(wb, tablename, overwrite = TRUE)
  
  
  
  # save figures into a pdf
  pdfname <- paste0("analysis_output/", nam, ".pdf")
  pdf(pdfname, height = 7, width = 11)
  for (plot in plot_list) {
    print(plot)
  }
  dev.off()
  
  
  # create height and width based on the number of clusters
  nclusters = length(unique(object@meta.data$cluster_id))
  
  # also plot individually for figures
  pdf(paste0("analysis_output/", nam ,"/umap.pdf"), height = 6, width = 6)
  print(umap_plot)
  dev.off()
  
  pdf(paste0("analysis_output/", nam ,"/umap2.pdf"), height = 6, width = 6)
  print(umap2)
  dev.off()
  
  pdf(paste0("analysis_output/", nam ,"/heat_map.pdf"), height = nclusters+3, width = nclusters+3)
  print(heat_map2)
  dev.off()
  
  pdf(paste0("analysis_output/", nam ,"/dot_plot.pdf"), height = nclusters/4+3, width = 5+length(specific_genes)/4)
  print(dot_plot2)
  dev.off()
  
  
  # also plot pcas of original and new objects
  p1 <- DimPlot(object, reduction = "pca", dims = c(1, 2), group.by = "ident") +
    ggtitle("PCA: PC1 vs PC2") + theme_minimal()
  
  p2 <- DimPlot(object, reduction = "pca", dims = c(1, 3), group.by = "ident") +
    ggtitle("PCA: PC1 vs PC3") + theme_minimal()
  
  p3 <- DimPlot(object, reduction = "pca", dims = c(2, 3), group.by = "ident") +
    ggtitle("PCA: PC2 vs PC3") + theme_minimal()
  
  p4 <- ElbowPlot(object) + 
    ggtitle("Elbow plot")
  
  pdf(paste0("analysis_output/sct/", nam ,"/", nam,"_pcas.pdf"), height = 7, width = 12)
  gridExtra::grid.arrange(umap_plot+ggtitle("original"), p1, NULL, p4, p2, p3, ncol = 3)
  
  
  p1 <- DimPlot(new_object, reduction = "pca", dims = c(1, 2), group.by = "ident") +
    ggtitle("PCA: PC1 vs PC2") + theme_minimal()
  
  p2 <- DimPlot(new_object, reduction = "pca", dims = c(1, 3), group.by = "ident") +
    ggtitle("PCA: PC1 vs PC3") + theme_minimal()
  
  p3 <- DimPlot(new_object, reduction = "pca", dims = c(2, 3), group.by = "ident") +
    ggtitle("PCA: PC2 vs PC3") + theme_minimal()
  
  p4 <- ElbowPlot(new_object) + 
    ggtitle("Elbow plot")
  
  gridExtra::grid.arrange(umap2+ggtitle("New"), p1, NULL, p4, p2, p3, ncol = 3)
  dev.off()
  
  pdf(paste0("analysis_output/", nam,"_color_testing.pdf"), height = 12, width = 15)
  gridExtra::grid.arrange(dot_plot1, dot_plot2, dot_plot3, nrow = 3)
  gridExtra::grid.arrange(heat_map1, heat_map2, heat_map3, ncol = 3)
  dev.off()
  
}






