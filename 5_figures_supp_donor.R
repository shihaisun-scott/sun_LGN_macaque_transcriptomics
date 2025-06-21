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



feature <- c(M = 2800, P = 3200, K = 1100)
resolutions <- c(M = 0.2, P = 0.26, K = 0.24)
ndims <- c(M = 4, P = 4, K = 7)


mpk_colors <- c("#de2d26", "#fee0d2",
                "#3182bd", "#a6bddb",
                "#2ca25f", "#99d8c9", "#e5f5f9")

exp_colors <- c("#2B7BB9", "#FF8C1A", "#3AAA35")
exp_colors2 <- c("#66c2a5", "#fc8d62", "#8da0cb")
exp_colors3 <- c("#1b9e77","#7570b3")
exp_colors4 <- c("#1b9e77","#d95f02","#7570b3")


# loop between each set of cells: plot UMAP but with species labeling
plot_list_donor <- list()
for (nam in names(subcells)) {
  object <- seurat_objects[[nam]]
  
  # plot UMAP
  dim_plot <- DimPlot(object, group.by = "donor", cols= exp_colors4, alpha = 0.85) + 
    ggtitle(nam) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme_void() 
  
  plot_list_donor[[nam]] <- dim_plot
  
}


plot_list_species <- list()

# loop between each set of cells: plot UMAP but with species labeling
for (nam in names(subcells)) {
  object <- seurat_objects[[nam]]
  
  # plot UMAP
  dim_plot <- DimPlot(object, group.by = "species", cols= exp_colors3, alpha = 0.85) + 
    ggtitle(nam) + 
    theme(aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme_void() 
  
  plot_list_species[[nam]] <- dim_plot
  
  
}

plots_id <- plot_grid(plotlist = c(plot_list_donor, plot_list_species), ncol = 3)

# save
pdfname <- paste0("analysis_output/5_supp_mpk_species_donor.pdf")
pdf(pdfname, height = 6, width = 9)
print(plots_id)
dev.off()



### plot K new vs K Bakken clustering
object <- seurat_objects[["K"]]
k_colors <- c("#2ca25f","#99d8c9","#e5f5f9")
object@meta.data$cluster_label_new <- newclusters$cluster_label[match(rownames(object@meta.data), newclusters$sample_name)]
object@meta.data$cluster_label_old <- metalist$cluster_label[match(rownames(object@meta.data), metalist$sample_name)]

# find k3 (new) vs K2 old
k_match = sum(object@meta.data$cluster_label_new == "K3" & object@meta.data$cluster_label_old == "K2")
k2_length = sum(object@meta.data$cluster_label_old == "K2")
perc_match = k_match / k2_length * 100
cat(perc_match, "%", "(", k_match, "/", k2_length, ")")

# new clustering
dim_plot_new <- DimPlot(object, group.by = "cluster_label_new", cols = c("#2B7BB9", "#FF8C1A", "#3AAA35"), alpha = 0.85) + 
  ggtitle(nam) + 
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5)
  )

# old clustering as second subplot
dim_plot_old <- DimPlot(object, group.by = "cluster_label_old", cols = c("#2B7BB9", "#3AAA35"), alpha = 0.85) + 
  ggtitle(nam) + 
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5)
  )

p_comb <- dim_plot_new + dim_plot_old
pdf("analysis_output/5_supp_k_new_old_compare.pdf", width = 6, height = 3)  # adjust size as needed
print(p_comb)
dev.off()


