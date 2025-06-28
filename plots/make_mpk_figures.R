# make_custom_figures.R
# note, journal full page width typically 7 inches

rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(cowplot)   
library(patchwork) 

theme(text = element_text(family = "Arial"))

fig_width <- c("M" = 7, "P" = 7, "K" = 7)
fig_height <- c("M" = 7, "P" = 7, "K" = 8)
volc_xlim <- list("M" = c(0, 2), "P" = c(0, 2), "K" = c(0, 10))

label_font_size = 12; # A B C ... font size

# load all of the saved plots
elbow_plot <- readRDS("elbow_plot.rds")
heat_map   <- readRDS("heat_map.rds")
dot_plot <- readRDS("dot_plot.rds")
dot_plot_vert <- readRDS("dot_plot_vert.rds")
gam_heat_map   <- readRDS("gam_heat_map.rds")
k_new_vs_bakken <- readRDS("k_new_vs_bakken.rds")
plot_list_donor   <- readRDS("plot_list_donor.rds")
plot_list_species <- readRDS("plot_list_species.rds")
umap_plot   <- readRDS("umap_plot.rds")
umap_plot_mpk   <- readRDS("umap_plot_mpk.rds")
umap_plot2   <- readRDS("umap_plot2.rds")
volcano_plot   <- readRDS("volcano_plot.rds")


# first set font sizes and remove legend if necessary
# define font sizes and remove title
custom_font <- function(p, base = 6) {
  p + theme(
    axis.title = element_text(size = base),
    axis.text = element_text(size = base - 1),
    legend.text = element_text(size = base - 1),
    legend.title = element_text(size = base),
    strip.text = element_text(size = base),
    plot.title = element_blank(),
    legend.position = "none"
  )
}

# x and y labels and ticks
remove_xy_all <- function(p) {
  p + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
}

# keep xy labels
# x and y labels and ticks
remove_xy_ticks <- function(p) {
  p + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
  )
}

# change the to a classic theme
elbow_plot <- lapply(elbow_plot, function(p) p + theme_bw())
elbow_plot <- lapply(elbow_plot, function(p) p + theme(axis.text.y = element_blank())) # remove y tick labels


# apply custom font and cleaning
elbow_plot <- lapply(elbow_plot, custom_font)
heat_map <- lapply(heat_map, custom_font)
dot_plot <- lapply(dot_plot, custom_font)
dot_plot_vert <- lapply(dot_plot_vert, custom_font)
gam_heat_map <- lapply(gam_heat_map, custom_font)
k_new_vs_bakken <- lapply(k_new_vs_bakken, custom_font)
plot_list_donor <- lapply(plot_list_donor, custom_font)
plot_list_species <- lapply(plot_list_species, custom_font)
umap_plot <- lapply(umap_plot, custom_font)
umap_plot_mpk <- lapply(umap_plot_mpk, custom_font)
umap_plot2 <- lapply(umap_plot2, custom_font)
volcano_plot <- lapply(volcano_plot, custom_font)

# remove xy labels and ticks
gam_heat_map <- lapply(gam_heat_map, remove_xy_ticks)
umap_plot <- lapply(umap_plot, remove_xy_all)
umap_plot2 <- lapply(umap_plot2, remove_xy_all)

## legend editing
gam_heat_map <- lapply(gam_heat_map, function(p) {
  p +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 5),      # Make text smaller
      legend.title = element_text(size = 6),            # Remove legend title
      legend.key.height = unit(0.3, "cm"),       # Smaller legend key height
      legend.key.width = unit(0.2, "cm")         # Smaller legend key width
    ) +
    scale_fill_viridis_c(
      labels = function(x) round(x, 2),
      name = "Sil.\nScore")  # Round labels to 1 decimal
})

# volcano legend on
volcano_plot <- lapply(volcano_plot, function(p) p + theme(legend.position = "right"))

# legend on right for dotplot and heatmap
legend_custom <- function(p) {
  p + theme(
    legend.position = "right",
    legend.text = element_text(size = 5),     
    # legend.title = element_text(size = 6),
    legend.title = element_blank(),
    # legend.title.position = "top",
    legend.key.height = unit(0.5, "cm"),       
    legend.key.width = unit(0.3, "cm")         
  )
}


dot_plot_vert <- lapply(dot_plot_vert, legend_custom)
heat_map <- lapply(heat_map, legend_custom)


# gam_heat_map <- lapply(gam_heat_map, function(p) p + theme(legend.position = "right"))

# square these
# elbow_plot <- lapply(elbow_plot, function(p) p + coord_fixed() + theme(aspect.ratio = 1))
# gam_heat_map <- lapply(gam_heat_map, function(p) p + coord_fixed() + theme(aspect.ratio = 1))


### VERSION 1- ALL
# the order is elbow, gam heat map, umap, umap2, heatmap, dotplot, volcano plot
### M PLOTTING
# Cell types to process
cell_types <- c("M", "P", "K")

# Create output directory
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# Loop through each cell type and build composite figure
for (ct in cell_types) {
  message("🔧 Building figure for: ", ct)
  
  # adjust volcano plot text size
  volcano_plot[[ct]]$layers[[4]]$aes_params$size <- 1.5
  volcano_plot[[ct]]$layers[[5]]$aes_params$size <- 1.5
  volcano_plot[[ct]] <- volcano_plot[[ct]] + xlim(volc_xlim[[ct]])
  
  elb_gam <- plot_grid(elbow_plot[[ct]],
                       gam_heat_map[[ct]],
                       nrow = 2,
                       rel_heights = c(1, 1))
  
  row1 <- plot_grid(
    elb_gam,
    plot_spacer(),
    umap_plot[[ct]],
    umap_plot2[[ct]],
    ncol = 4,
    rel_widths = c(0.9, 0.2, 1, 1),
    labels = c("A", "", "B", "C"),
    label_size = label_font_size,
    label_fontface = "bold"
  )
  
  row2 <- plot_grid(
    heat_map[[ct]],
    dot_plot_vert[[ct]],
    ncol = 2,
    rel_widths = c(1.2, 0.8),
    labels = c("D", "E"),
    label_size = label_font_size,
    label_fontface = "bold"
  )
  
  row3 <- plot_grid(
    volcano_plot[[ct]],
    ncol = 1,
    rel_widths = c(0.8),
    labels = "F",
    label_size = label_font_size
  )
  
  fig_combined <- plot_grid(
    row1, row2, row3,
    nrow = 3,
    rel_heights = c(1, 1.2, 0.8)
  )
  
  save_name <- paste0("Figures/Figure_", ct)
  ggsave(
    filename = paste0(save_name, ".pdf"),
    plot = fig_combined,
    width = fig_width[[ct]],
    height = fig_height[[ct]],
    units = "in"
  )
  
  ggsave(
    filename = paste0(save_name, ".png"),
    plot = fig_combined,
    width = fig_width[[ct]],
    height = fig_height[[ct]],
    units = "in",
    dpi = 600   
  )
  
  # ggsave(
  #   filename = paste0(save_name, ".svg"),
  #   plot = fig_combined,
  #   width = fig_width[[ct]],
  #   height = fig_height[[ct]],
  #   units = "in"
  # )
}



