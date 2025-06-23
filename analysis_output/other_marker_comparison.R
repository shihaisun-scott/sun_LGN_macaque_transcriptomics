# compare_markers.R

# Load libraries
library(dplyr)
library(openxlsx)
library(tools)

# --- Paths ---

# Folder with new marker files
new_marker_path <- "analysis_output"

# Folder with old marker files
old_marker_path <- "data/first_version_stats"

# Output file
output_file <- "analysis_output/gene_comparison_summary.xlsx"

# --- Load marker tables ---

load_marker_tables <- function(folder_path) {
  files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)
  marker_list <- list()
  
  for (f in files) {
    wb <- loadWorkbook(f)
    if (!"Summary" %in% names(wb)) next  # Skip files with no Summary sheet
    marker_data <- read.xlsx(f, sheet = "Summary")
    
    if (!("top_genes" %in% colnames(marker_data)) || !("cluster" %in% colnames(marker_data))) next
    
    for (i in 1:nrow(marker_data)) {
      cluster <- as.character(marker_data$cluster[i])
      genes <- strsplit(marker_data$top_genes[i], ",\\s*")[[1]]
      marker_list[[cluster]] <- data.frame(
        gene = genes,
        cluster = cluster,
        stringsAsFactors = FALSE
      )
    }
  }
  return(marker_list)
}

# --- Load new and old markers ---

new_markers_list <- load_marker_tables(new_marker_path)
old_marker_files <- list.files(old_marker_path, pattern = "^Cluster_.*\\.xlsx$", full.names = TRUE)

old_marker_list <- list()
for (f in old_marker_files) {
  cluster <- gsub("Cluster_|\\.xlsx", "", basename(f))
  old_df <- read.xlsx(f)
  if (all(c("gene", "p_val_adj", "avg_log2FC") %in% colnames(old_df))) {
    old_marker_list[[cluster]] <- old_df %>%
      select(gene, p_val_adj, avg_log2FC) %>%
      rename(old_p_val_adj = p_val_adj, old_log2FC = avg_log2FC)
  }
}

# --- Compare markers ---

comparison_all <- list()

for (cl in names(new_markers_list)) {
  new_genes <- new_markers_list[[cl]] %>%
    mutate(new_p_val_adj = NA, new_log2FC = NA)
  
  # Try to load detailed markers (if available)
  new_marker_file <- file.path(new_marker_path, paste0("Cluster_", cl, ".xlsx"))
  if (file.exists(new_marker_file)) {
    detailed <- read.xlsx(new_marker_file)
    if (all(c("gene", "p_val_adj", "avg_log2FC") %in% colnames(detailed))) {
      new_genes <- left_join(new_genes, detailed %>%
                               select(gene, p_val_adj, avg_log2FC) %>%
                               rename(new_p_val_adj = p_val_adj, new_log2FC = avg_log2FC),
                             by = "gene")
    }
  }
  
  old_genes <- old_marker_list[[cl]]
  
  comparison <- full_join(new_genes, old_genes, by = "gene")
  
  comparison <- comparison %>%
    mutate(status = case_when(
      !is.na(new_log2FC) & is.na(old_log2FC) ~ "New",
      !is.na(new_log2FC) & !is.na(old_log2FC) ~ "Common",
      is.na(new_log2FC) & !is.na(old_log2FC) ~ "Missing",
      TRUE ~ "NA"
    )) %>%
    mutate(cluster = cl) %>%
    select(cluster, gene, status, new_log2FC, new_p_val_adj, old_log2FC, old_p_val_adj)
  
  comparison_all[[cl]] <- comparison
}

# Combine all clusters
comparison_df <- bind_rows(comparison_all)

# --- Save output ---

wb <- createWorkbook()
addWorksheet(wb, "Gene_Comparison")
writeData(wb, "Gene_Comparison", comparison_df)
saveWorkbook(wb, output_file, overwrite = TRUE)

message("Comparison saved to: ", output_file)
