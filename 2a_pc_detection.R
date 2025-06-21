# Elbow plots for M, P, K — variance explained for increasing feature steps
# Red dotted lines placed on elbow (qualitative assessment)

rm(list = ls())


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(Matrix)

# load data from 1_data_collection
load("data/data.rda")
load("data/precluster_info.rda")


subdat=datlist[,precluster$sample_name]

# Categorize cells into MPK
plot_objects=list()
subcells=list()
subcells[["M"]] <- subset(metalist$sample_name, precluster$cluster_label == "M")
subcells[["P"]] <- subset(metalist$sample_name, precluster$cluster_label == "P")
subcells[["K"]] <- subset(metalist$sample_name, grepl("^K", precluster$cluster_label))

# Define range of feature steps to test
feature_steps <- seq(1000, 4000, by = 500)

# Store variance explained per group
pc_variances_all <- list()

for (nam in names(subcells)) {
  cat("Processing group:", nam, "\n")
  
  keepcells <- subcells[[nam]]
  expr <- subdat[, keepcells]
  cell_meta <- newclusters[keepcells, , drop = FALSE]
  
  pc_variances <- list()
  
  for (f in feature_steps) {
    cat("  Features:", f, "\n")
    
    # Create object
    obj_tmp <- CreateSeuratObject(expr, meta.data = cell_meta)
    obj_tmp[["RNA"]] <- split(obj_tmp[["RNA"]], f = obj_tmp$donor)
    obj_tmp <- SCTransform(object = obj_tmp, verbose = FALSE, return.only.var.genes = FALSE,
                          variable.features.n = f)
    obj_tmp <- RunPCA(obj_tmp, npcs = 30, verbose = FALSE)
    
    # get variance per PC
    sdev <- obj_tmp[["pca"]]@stdev
    var_explained <- sdev^2 / sum(sdev^2)
    pc_variances[[as.character(f)]] <- var_explained
  }
  
  # assemble into dataframe
  pc_df <- do.call(rbind, lapply(names(pc_variances), function(f) {
    data.frame(
      PC = 1:length(pc_variances[[f]]),
      Variance = pc_variances[[f]],
      Features = as.numeric(f),
      Group = nam
    )
  }))
  
  pc_variances_all[[nam]] <- pc_df
}


## plot with overlapping lines
# combine
combined_pc_df <- bind_rows(pc_variances_all)
combined_pc_df$Features <- factor(combined_pc_df$Features)

# Plot elbow curves: one plot per group, overlapping feature lines
elb_plots <- ggplot(combined_pc_df, aes(x = PC, y = Variance, group = Features, color = Features)) +
  geom_line(alpha = 0.5) +
  geom_vline(xintercept = c(4,4,7), linetype = "dashed", color = "red") +
  facet_wrap(~ Group, ncol = 3) +
  labs(title = "Elbow Plots for M, P, K - Varying Feature Counts",
       x = "Principal Component",
       y = "Proportion of Variance Explained") +
  theme_minimal()

# Save to PDF
pdfname <- "2_pc_detection_overlap.pdf"
pdf(pdfname, height = 5, width = 14)
print(elb_plots)
dev.off()

# Show the plot
elb_plots


