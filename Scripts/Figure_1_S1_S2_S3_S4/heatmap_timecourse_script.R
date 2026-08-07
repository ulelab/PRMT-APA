library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)

# =============================================================================
# Setup — read files and compute metrics

###setwd
setwd("Data/Figure_1_S2_S3_S4/APA_CSV_files/timecourse/positional_info")

csv_files     <- list.files(pattern = "\\.csv$")
all_sites_csv <- lapply(csv_files, read_csv)

# Add proximal-to-distal shift column to every table
all_sites_csv <- lapply(all_sites_csv, function(df) {
  df %>%
    group_by(gene_name) %>%
    mutate(
      prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] -
        change_in_usage[pA_type == "proximal pA"]
    ) %>%
    ungroup()
})

# Split into sig (1:7) and all-sites (8:14) lists
sig_sites_list <- all_sites_csv[1:7]
all_sites_list <- all_sites_csv[8:14]

# =============================================================================
# Step 1: Build the "sig anchor" feature set
# For each condition, keep distal pA rows with |prox_to_dist_shift| > 0.3,
# then take the union of sig features across all conditions.
# This replaces the 7 repetitive anchor-then-merge chains.
# =============================================================================

join_keys <- c("gene_name", "feature_id", "chr", "start", "end", "strand")

sig_features <- lapply(sig_sites_list, function(df) {
  df %>%
    filter(abs(prox_to_dist_shift) > 0.3, pA_type == "distal pA") %>%
    select(all_of(c(join_keys, "UTR_type"))) %>%
    dplyr::slice(1, .by = gene_name)   # one row per gene
}) %>%
  bind_rows() %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  distinct(gene_name,  .keep_all = TRUE)

# =============================================================================
# Step 2: For each condition's all-sites table, extract change_in_usage
# (renamed with condition suffix) and join onto the sig feature set.
# =============================================================================

# Rename change_in_usage with condition suffix, keep only join keys + that column
prep_condition <- function(df) {
  condition_label <- df$condition[1]
  new_col         <- paste0(condition_label, "_change_in_usage")
  
  df %>%
    filter(pA_type == "distal pA") %>%
    dplyr::slice(1, .by = gene_name) %>%
    select(all_of(join_keys), change_in_usage) %>%
    dplyr::rename(!!new_col := change_in_usage)
}

condition_tables <- lapply(all_sites_list, prep_condition)

# Left-join each condition table onto the sig feature set in one chain
all_sig_sites <- condition_tables %>%
  purrr::reduce(function(acc, df) left_join(acc, df, by = join_keys),
         .init = sig_features)

# =============================================================================
# Step 3: Filter and prepare for heatmap
# (same logic as original)
# =============================================================================

# Drop any gene with an NA in any condition column
all_sig_sites <- all_sig_sites %>%
  filter(rowSums(is.na(.)) == 0)

# Split by UTR type
TUTR_sig_sites     <- all_sig_sites %>% filter(UTR_type == "TUTR")
non_TUTR_sig_sites <- all_sig_sites %>% filter(!UTR_type %in% c("TUTR", "iAPA"))

# Extract change_in_usage columns and convert to % for heatmap
extract_heatmap_mat <- function(df) {
  df %>%
    dplyr::select(ends_with("_change_in_usage")) %>%
    mutate(across(everything(), ~ .x * 100))
}

TUTR_heatmap     <- extract_heatmap_mat(TUTR_sig_sites)
non_TUTR_heatmap <- extract_heatmap_mat(non_TUTR_sig_sites)

# =============================================================================
# Step 4: Plot
# =============================================================================

RdBu   <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(40)
breaks <- seq(-40, 40, length.out = 41)

TUTR_heatmap <- TUTR_heatmap %>%
  dplyr::select(
    DMAi_24hrs_change_in_usage,
    DMAi_48hrs_change_in_usage,
    DMAi_72hrs_change_in_usage,
    DMAi_96hrs_change_in_usage,
    Palbo_24hr_change_in_usage,
    Palbo_48hr_change_in_usage,
    CDK11i_24hrs_change_in_usage
  )

non_TUTR_heatmap <- non_TUTR_heatmap %>%
  dplyr::select(
    DMAi_24hrs_change_in_usage,
    DMAi_48hrs_change_in_usage,
    DMAi_72hrs_change_in_usage,
    DMAi_96hrs_change_in_usage,
    Palbo_24hr_change_in_usage,
    Palbo_48hr_change_in_usage,
    CDK11i_24hrs_change_in_usage
  )

col_labels <- c(
  "DMAi_24hrs_change_in_usage" = "DMAi 24hrs",
  "DMAi_48hrs_change_in_usage" = "DMAi 48hrs",
  "DMAi_72hrs_change_in_usage" = "DMAi 72hrs",
  "DMAi_96hrs_change_in_usage" = "DMAi 96hrs",
  "CDK11i_24hrs_change_in_usage" = "CDK11i/TOPKi 24hrs",
  "Palbo_24hr_change_in_usage" = "CDK4i/CDK6i 24hrs",
  "Palbo_48hr_change_in_usage" = "CDK4i/CDK6i 48hrs"
)

colnames(TUTR_heatmap) <- col_labels[colnames(TUTR_heatmap)]
colnames(non_TUTR_heatmap) <- col_labels[colnames(non_TUTR_heatmap)]

pheatmap(
  non_TUTR_heatmap,
  clustering_distance_rows = "correlation",
  cluster_cols = FALSE,
  color                    = RdBu,
  breaks                   = breaks,
  show_rownames            = FALSE,
  fontsize_col             = 14,
  angle_col                = "45",
  legend_breaks            = c(-50, -25, 0, 25, 50),
  legend_labels            = c("-50%", "-25%", "0%", "25%", "50%"),
  cutree_rows              = 2
)
