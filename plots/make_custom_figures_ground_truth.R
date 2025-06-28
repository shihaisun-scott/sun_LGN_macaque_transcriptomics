# make custom figure for the ground truth data:
# 2x2 plot
# initial two (left), then clustering two (right)

rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(cowplot)   
library(patchwork) 

theme(text = element_text(family = "Arial"))

label_font_size = 12; # A B C ... font size

# load plot objects
cluster_umap1 <- readRDS("GT_1_cluster.rds")
truth_umap1 <- readRDS("GT_1_truth.rds")
cluster_umap2 <- readRDS("GT_2_cluster.rds")
truth_umap2 <- readRDS("GT_2_truth.rds")

# define font sizes 
custom_font <- function(p, base = 6) {
  p + theme(
    axis.title = element_text(size = base),
    axis.text = element_text(size = base - 1),
    legend.text = element_text(size = base - 1),
    legend.title = element_text(size = base),
    strip.text = element_text(size = base),
    plot.title = element_blank()
  )
}

cluster_umap1 <- custom_font(cluster_umap1)
truth_umap1 <- custom_font(truth_umap1)
cluster_umap2 <- custom_font(cluster_umap2)
truth_umap2 <- custom_font(truth_umap2)

# adjust dot sizes
adjust_dot_size <- function(p, size = 0.5) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
      p$layers[[i]]$aes_params$size <- size
    }
  }
  return(p)
}
cluster_umap1 <- adjust_dot_size(cluster_umap1, size = 0.05)
truth_umap1 <- adjust_dot_size(truth_umap1, size = 0.05)
cluster_umap2 <- adjust_dot_size(cluster_umap2, size = 0.5)
truth_umap2 <- adjust_dot_size(truth_umap2, size = 0.5)

row1 <- plot_grid(
  cluster_umap1,
  cluster_umap2,
  ncol = 2,
  rel_widths = c(1,1),
  labels = c("A", "B"),
  label_size = label_font_size,
  label_fontface = "bold"
)

row2 <- plot_grid(
  truth_umap1,
  truth_umap2,
  ncol = 2,
  rel_widths = c(1, 1),
  # labels = c("", ""),
  label_size = label_font_size,
  label_fontface = "bold"
)

fig_combined <- plot_grid(
  row1, row2,
  nrow = 2,
  rel_heights = c(1,1)
)

ggsave(
  filename = "Figures/Supp_Figure_GT.pdf",
  plot = fig_combined,
  width = 9,
  height = 5,
  units = "in"
)

ggsave(
  filename = "Figures/Supp_Figure_GT.png",
  plot = fig_combined,
  width = 9,
  height = 5,
  units = "in",
  dpi = 600
)



