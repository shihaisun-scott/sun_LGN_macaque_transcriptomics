# make two supp figures:
# 1. donor and species for mpk
# 2. old and new labels for k clusters


rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(cowplot)   
library(patchwork) 

theme(text = element_text(family = "Arial"))

label_font_size = 12; # A B C ... font size

# load plot objects
plot_list_donor <- readRDS("plot_list_donor.rds")
plot_list_species <- readRDS("plot_list_species.rds")
k_new_vs_bakken <- readRDS("k_new_vs_bakken.rds")

# define font sizes 
custom_font <- function(p, base = 6) {
  p + theme(
    axis.title = element_text(size = base),
    axis.text = element_text(size = base - 1),
    legend.text = element_text(size = base - 1),
    legend.title = element_text(size = base),
    strip.text = element_text(size = base),
    plot.title = element_text(size = base + 2, face = "plain")
  )
}

plot_list_donor <- lapply(plot_list_donor, custom_font)
plot_list_species <- lapply(plot_list_species, custom_font)
k_new_vs_bakken <- lapply(k_new_vs_bakken, custom_font)

# adjust dot sizes
adjust_dot_size <- function(p, size = 0.5) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
      p$layers[[i]]$aes_params$size <- size
    }
  }
  return(p)
}

plot_list_donor <- lapply(plot_list_donor, adjust_dot_size, size = 0.5)
plot_list_species <- lapply(plot_list_species, adjust_dot_size, size = 0.5)
k_new_vs_bakken <- lapply(k_new_vs_bakken, adjust_dot_size, size = 0.5)


# donor and species first
row1 <- plot_grid(
  plotlist = plot_list_donor,
  ncol = 3,
  rel_widths = c(1,1,1)
  
)

row2 <- plot_grid(
  plotlist = plot_list_species,
  ncol = 3,
  rel_widths = c(1,1,1)
)

fig_combined <- plot_grid(
  row1, row2,
  nrow = 2,
  rel_heights = c(1,1),
  labels = "AUTO",
  label_size = label_font_size,
  label_fontface = "bold"
)

ggsave(
  filename = "Figures/Supp_Figure_species.pdf",
  plot = fig_combined,
  width = 7,
  height = 3,
  units = "in"
)

ggsave(
  filename = "Figures/Supp_Figure_species.png",
  plot = fig_combined,
  width = 7,
  height = 3,
  units = "in",
  dpi = 600
)


## do the k plots
fig_combined <- plot_grid(
  plot_list = k_new_vs_bakken,
  ncol = 2,
  labels = "AUTO",
  label_size = label_font_size,
  label_fontface = "bold"
)

ggsave(
  filename = "Figures/Supp_Figure_K_bakken_compare.pdf",
  plot = fig_combined,
  width = 7,
  height = 3,
  units = "in"
)

ggsave(
  filename = "Figures/Supp_Figure_K_bakken_compare.png",
  plot = fig_combined,
  width = 7,
  height = 3,
  units = "in",
  dpi = 600
)

