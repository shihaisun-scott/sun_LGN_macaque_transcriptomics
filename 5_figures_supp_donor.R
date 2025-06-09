# supp 1: plot each of MPK with animal id labeled
# supp 2: plot K with new clustering versus old clustering (Kap and Kp)

rm(list = ls())


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(patchwork)

load("data/data.rda")
load("data/precluster_info.rda")
load("data/sct_newclusters.rda")
load("data/sct_seurat_objects.rda")

subdat=datlist[,precluster$sample_name]
plot_objects=list()
subcells=list()
subcells[["M"]] <- subset(metalist$sample_name, precluster$cluster_label == "M")
subcells[["P"]] <- subset(metalist$sample_name, precluster$cluster_label == "P")
subcells[["K"]] <- subset(metalist$sample_name, grepl("^K", precluster$cluster_label))
# subcells[["K1"]] <- subset(metalist$sample_name, precluster$cluster_label == "K1")
# subcells[["K2"]] <- subset(metalist$sample_name, precluster$cluster_label == "K2")



# custom settings
feature <- c(M = 30, P = 40, K = 220)
ndims <- c(M = 3, P = 3, K = 3)
resolutions <- c(M = 0.32, P = 0.2, K = 0.28)
exp_colors <- c("purple", "black", "yellow")

plot_list <- list()

# loop between each set of cells: plot UMAP but with species labeling
for (nam in names(subcells)) {
  object <- seurat_objects[[nam]]
  
  # plot UMAP
  dim_plot <- DimPlot(object, group.by = "donor", cols= exp_colors, alpha = 0.5) + 
    ggtitle(nam) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme_void() 
  
  plot_list[[nam]] <- dim_plot
  
  
}

# save
pdfname <- paste0("analysis_output/supp_mpk_donor_id.pdf")
pdf(pdfname, height = 6, width = 6)
for (plot in plot_list) {
  
  print(plot)
}
dev.off()

plot_list <- list()

# loop between each set of cells: plot UMAP but with species labeling
for (nam in names(subcells)) {
  object <- seurat_objects[[nam]]
  
  # plot UMAP
  dim_plot <- DimPlot(object, group.by = "donor", cols= exp_colors, alpha = 0.5) + 
    ggtitle(nam) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme_void() 
  
  plot_list[[nam]] <- dim_plot
  
  
}

# save
pdfname <- paste0("analysis_output/supp_mpk_donor_id.pdf")
pdf(pdfname, height = 6, width = 6)
for (plot in plot_list) {
  
  print(plot)
}
dev.off()


plot_list_species <- list()

# loop between each set of cells: plot UMAP but with species labeling
for (nam in names(subcells)) {
  object <- seurat_objects[[nam]]
  
  # plot UMAP
  dim_plot <- DimPlot(object, group.by = "species", cols= exp_colors, alpha = 0.5) + 
    ggtitle(nam) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme_void() 
  
  plot_list_species[[nam]] <- dim_plot
  
  
}

# save
pdfname <- paste0("analysis_output/supp_mpk_species_id.pdf")
pdf(pdfname, height = 6, width = 6)
for (plot in plot_list_species) {
  
  print(plot)
}
dev.off()



### plot K new vs K Bakken clustering
object <- seurat_objects[["K"]]
k_colors <- colorRampPalette(colors = c("#004B49","#ACE1AF"))(3)
object@meta.data$cluster_label_new <- newclusters$cluster_label[match(rownames(object@meta.data), newclusters$sample_name)]
object@meta.data$cluster_label_old <- metalist$cluster_label[match(rownames(object@meta.data), metalist$sample_name)]

# new clustering
dim_plot_new <- DimPlot(object, group.by = "cluster_label_new", cols = k_colors, alpha = 0.5) + 
  ggtitle(nam) + 
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5)
  )

# old clustering as second subplot
dim_plot_old <- DimPlot(object, group.by = "cluster_label_old", cols = c("#004B49","#ACE1AF"), alpha = 0.5) + 
  ggtitle(nam) + 
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5)
  )

p_comb <- dim_plot_new + dim_plot_old
pdf("analysis_output/supp_k_new_old_compare.pdf", width = 10, height = 5)  # adjust size as needed
print(p_comb)
dev.off()


