### Clear environment
rm(list = ls())

### Load libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(gridExtra)
library(cluster)
library(mgcv)
library(clue)

### Load clustered Seurat object
obj <- readRDS("pbmc_clustered.rds")

### Subset cluster 0
cluster0_cells <- WhichCells(obj, idents = "0")
cluster0_cells <- cluster0_cells[obj@meta.data[cluster0_cells, "cell_label"] != "cytotoxic_t"]
cluster0_cells <- cluster0_cells[obj@meta.data[cluster0_cells, "cell_label"] != "cd14_monocytes"]

expr <- GetAssayData(obj[, cluster0_cells], slot = "counts")
cell_labels_cluster0 <- obj$cell_label[cluster0_cells]

### Set parameters
# sequential optimisation with 87% match
# nfeatures = 220
# ndim = 6
# res = 0.43

# brute optim with 91% match
# nfeatures = 3000
# ndim = 3
# res = 0.48

nfeatures = 3000
ndim = 7
res = 0.7


# Number of iterations
n_iterations <- 100

# Store match percentages
match_percentages <- numeric(n_iterations)

# Begin bootstrapping
prev_cells <- NULL
for (i in 1:n_iterations) {
  cat(sprintf("Bootstrap iteration: %d\n", i))
  set.seed(123 + i)
  # Randomly remove 20% of cells
  sampled_cells <- sample(cluster0_cells, size = floor(0.8 * length(cluster0_cells)))
  
  # check if same cells
  if (!is.null(prev_cells)) {
    changed <- sum(!(sampled_cells %in% prev_cells))
    cat(sprintf("Changed cells from previous: %d\n", changed))
  }
  prev_cells <- sampled_cells
  
  expr_sub <- GetAssayData(obj[, sampled_cells], slot = "counts")
  labels_sub <- obj$cell_label[sampled_cells]
  
  # Seurat analysis pipeline
  object <- CreateSeuratObject(expr_sub)
  object$cell_label <- labels_sub
  object <- SCTransform(object, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = nfeatures)
  object <- RunPCA(object, npcs = 20, verbose = FALSE)
  object <- RunUMAP(object, dims = 1:ndim, verbose = FALSE)
  object <- FindNeighbors(object, dims = 1:ndim, verbose = FALSE)
  object <- FindClusters(object, resolution = res, verbose = FALSE)
  
  # Calculate match percentage
  # clusters <- Idents(object)
  # labels <- object$cell_label
  # confusion <- table(clusters, labels)
  # correct_matches <- sum(apply(confusion, 1, max))
  # total_cells <- length(clusters)
  # percent_match <- correct_matches / total_cells * 100
  # match_percentages[i] <- percent_match
  # 
  clusters <- Idents(object)
  labels <- object$cell_label
  confusion <- table(clusters, labels)
  conf_mat <- as.matrix(confusion)
  n <- max(nrow(conf_mat), ncol(conf_mat))
  square_mat <- matrix(0, n, n)
  square_mat[1:nrow(conf_mat), 1:ncol(conf_mat)] <- conf_mat
  cost_mat <- max(square_mat) - square_mat
  assignment <- solve_LSAP(cost_mat)
  correct_matches <- sum(square_mat[cbind(1:n, assignment)])
  total_cells <- length(clusters)
  percent_match <- correct_matches / total_cells * 100
  cat(sprintf(">> Optimally matched percentage: %.2f%%\n", percent_match))
  match_percentages[i] <- percent_match
}

# Final output
mean_match <- mean(match_percentages)
sd_match <- sd(match_percentages)
cat(sprintf(">> Mean match over %d iterations: %.2f%% (SD = %.2f%%)\n", n_iterations, mean_match, sd_match))


