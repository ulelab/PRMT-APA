# =============================================================================
# PCF11 overexpression rescue of DMAi-induced 3'UTR lengthening (LU-99)
# TUTR genes only. Clones c3 (30ng dox) and c7 (20ng dox); GFP (5ng) control.
# =============================================================================

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(ggpattern)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)


# =============================================================================
# SECTION 1: LOAD AND MERGE
# =============================================================================

PCF11_csv_dir <- "/Users/k2362866/Documents/AZ_postdoc/PRMT-APA_resubmission/Data/Figure_4_S7/APA_CSV_files/APA_classified/no_CSTF2"

csv_files <- list.files(
  path = PCF11_csv_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

PCF11_oe_files <- lapply(csv_files, function(file) {
  df <- read.csv(file)
  
  df <- df %>%
    dplyr::group_by(gene_name) %>%
    dplyr::mutate(
      prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] -
        change_in_usage[pA_type == "proximal pA"]
    ) %>%
    dplyr::ungroup()
  
  df <- df %>%
    dplyr::filter(UTR_type == "TUTR") %>%
    dplyr::select(
      gene_name, feature_id, chr, start, end, strand, UTR_type,
      pA_usage_control, pA_usage_treatment, change_in_usage,
      prox_to_dist_shift, twostep_feature_padj, condition,
      pA_type, mean_usage
    )
  
  condition_value <- unique(df$condition)[1]
  
  columns_to_rename <- c(
    "prox_to_dist_shift", "twostep_feature_padj", "change_in_usage",
    "mean_usage", "pA_usage_control", "pA_usage_treatment"
  )
  
  df %>%
    dplyr::rename_with(
      .cols = all_of(columns_to_rename),
      .fn = ~ paste0(condition_value, "_", .x)
    )
})

merge_dfs <- function(df1, df2) {
  dplyr::inner_join(
    df1, df2,
    by = c("gene_name", "feature_id", "chr", "start", "end",
           "strand", "pA_type", "UTR_type")
  )
}

merged_df <- Reduce(merge_dfs, PCF11_oe_files)

# keep only genes with a clean proximal/distal pair
merged_df <- merged_df %>%
  dplyr::group_by(gene_name) %>%
  dplyr::filter(
    dplyr::n() == 2,
    any(pA_type == "proximal pA"),
    any(pA_type == "distal pA")
  ) %>%
  dplyr::ungroup()

# =============================================================================
# SECTION 2: CONTROL PAIRS AND DMAi-LENGTHENED GENE SET
# =============================================================================

clone_config <- list(
  c3 = list(clone = "c3", dose = "30ng"),
  c7 = list(clone = "c7", dose = "20ng")
)

# DMAi-lengthened: shift > 0.2 in BOTH clones
dmai_definition <- merged_df %>%
  dplyr::distinct(
    gene_name,
    LU99_PCF11_c3_no_dox_DMSO_vs_no_dox_DMAi_prox_to_dist_shift,
    LU99_PCF11_c3_no_dox_DMSO_vs_no_dox_DMAi_twostep_feature_padj,
    LU99_PCF11_c7_no_dox_DMSO_vs_no_dox_DMAi_prox_to_dist_shift,
    LU99_PCF11_c7_no_dox_DMSO_vs_no_dox_DMAi_twostep_feature_padj
  )

dmai_lengthened <- dmai_definition %>%
  dplyr::filter(
    LU99_PCF11_c3_no_dox_DMSO_vs_no_dox_DMAi_prox_to_dist_shift > 0.2,
    LU99_PCF11_c7_no_dox_DMSO_vs_no_dox_DMAi_prox_to_dist_shift > 0.2,
    (
      LU99_PCF11_c3_no_dox_DMSO_vs_no_dox_DMAi_twostep_feature_padj < 0.05 |
        LU99_PCF11_c7_no_dox_DMSO_vs_no_dox_DMAi_twostep_feature_padj < 0.05
    )
  ) %>%
  dplyr::pull(gene_name)

# =============================================================================
# SECTION 3: MITIGATION RATIOS PER CLONE
# =============================================================================

compute_ratio <- function(merged_df, line, dose, keep_genes) {
  alone_col <- glue::glue("LU99_{line}_no_dox_DMSO_vs_no_dox_DMAi_prox_to_dist_shift")
  oe_col    <- glue::glue("LU99_{line}_no_dox_DMSO_vs_{dose}_dox_DMAi_prox_to_dist_shift")
  
  merged_df %>%
    dplyr::filter(gene_name %in% keep_genes) %>%
    dplyr::mutate(
      shift_alone = .data[[alone_col]],
      shift_oe    = .data[[oe_col]],
      ratio       = .data[[oe_col]] / .data[[alone_col]]
    )
}

c3_compare <- compute_ratio(merged_df, line = "PCF11_c3", dose = "30ng",
                            keep_genes = dmai_lengthened) %>%
  dplyr::rename(shift_alone_c3 = shift_alone,
                shift_oe_c3    = shift_oe,
                ratio_c3       = ratio)

c7_compare <- compute_ratio(merged_df, line = "PCF11_c7", dose = "20ng",
                            keep_genes = dmai_lengthened) %>%
  dplyr::rename(shift_alone_c7 = shift_alone,
                shift_oe_c7    = shift_oe,
                ratio_c7       = ratio)

gfp_compare <- compute_ratio(merged_df, line = "GFP", dose = "5ng",
                             keep_genes = dmai_lengthened) %>%
  dplyr::rename(shift_alone_gfp = shift_alone,
                shift_oe_gfp    = shift_oe,
                ratio_gfp       = ratio)

# =============================================================================
# SECTION 4: COMBINED TABLE (two rows per gene: proximal + distal)
# =============================================================================

usage_c3 <- "LU99_PCF11_c3_no_dox_DMSO_vs_no_dox_DMAi_pA_usage_treatment"
usage_c7 <- "LU99_PCF11_c7_no_dox_DMSO_vs_no_dox_DMAi_pA_usage_treatment"

combined <- c3_compare %>%
  dplyr::select(
    gene_name, feature_id, pA_type,
    chr, start, end, strand,
    shift_alone_c3, shift_oe_c3, ratio_c3,
    dmai_prox_usage_c3 = all_of(usage_c3)
  ) %>%
  dplyr::inner_join(
    c7_compare %>%
      dplyr::select(gene_name, feature_id,
                    shift_alone_c7, shift_oe_c7, ratio_c7,
                    dmai_prox_usage_c7 = all_of(usage_c7)),
    by = c("gene_name", "feature_id")
  ) %>%
  dplyr::inner_join(
    gfp_compare %>% dplyr::select(gene_name, feature_id, ratio_gfp, shift_oe_gfp),
    by = c("gene_name", "feature_id")
  ) %>%
  dplyr::mutate(
    dmai_prox_usage  = rowMeans(cbind(dmai_prox_usage_c3, dmai_prox_usage_c7), na.rm = TRUE),
    ratio_mean       = rowMeans(cbind(ratio_c3, ratio_c7), na.rm = TRUE),
    shift_alone_mean = rowMeans(cbind(shift_alone_c3, shift_alone_c7), na.rm = TRUE),
    shift_oe_mean    = rowMeans(cbind(shift_oe_c3, shift_oe_c7), na.rm = TRUE)
  )

# contiguous bins on the clone-mean ratio
combined <- combined %>%
  dplyr::mutate(
    mitigation_category = dplyr::case_when(
      ratio_mean <  0     ~ "reversed",
      ratio_mean <  0.25  ~ "full mitigation",
      ratio_mean <  0.6   ~ "partial mitigation",
      ratio_mean <  0.85  ~ "weak mitigation",
      ratio_mean <= 1.2   ~ "no mitigation",
      ratio_mean >  1.2   ~ "exacerbated",
      TRUE                ~ "undetermined"      # only catches NAs
    )
  )

cat_levels <- c("exacerbated", "no mitigation", "weak mitigation",
                "partial mitigation", "full mitigation", "reversed")

# ---- SOURCE OF TRUTH #1: two rows per gene (site-level) ----
# used for BED/FASTA, proximal-usage and distance plots
TUTR_categorized_both_rows <- combined %>%
  dplyr::filter(!is.na(mitigation_category), mitigation_category != "undetermined") %>%
  dplyr::mutate(
    mitigation_category = factor(mitigation_category, levels = cat_levels),
    motif_group = dplyr::case_when(
      mitigation_category %in% c("full mitigation", "reversed")                     ~ "strong_mitigation",
      mitigation_category %in% c("weak mitigation", "no mitigation", "exacerbated") ~ "not_mitigated",
      TRUE                                                                          ~ "other"  # partial
    )
  )

# ---- SOURCE OF TRUTH #2: one row per gene (gene-level) ----
# used for the heatmap only (prox_to_dist_shift is identical on both rows)
TUTR_data_categorized_single <- TUTR_categorized_both_rows %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(mitigation_category, ratio_mean)

# sanity check: same genes, half the rows
stopifnot(
  setequal(TUTR_data_categorized_single$gene_name,
           unique(TUTR_categorized_both_rows$gene_name))
)

# =============================================================================
# SECTION 5: HEATMAPS (gene-level)
# =============================================================================

heatmap_data <- as.matrix(
  TUTR_data_categorized_single[, c("shift_alone_mean", "shift_oe_c3", "shift_oe_c7")]
)

heatmap_data_with_gfp <- as.matrix(
  TUTR_data_categorized_single[, c("shift_alone_mean", "shift_oe_c3",
                                   "shift_oe_c7", "shift_oe_gfp")]
)

annotations <- data.frame(
  Category = TUTR_data_categorized_single$mitigation_category
)
rownames(annotations)           <- TUTR_data_categorized_single$gene_name
rownames(heatmap_data)          <- TUTR_data_categorized_single$gene_name
rownames(heatmap_data_with_gfp) <- TUTR_data_categorized_single$gene_name

annotation_colors <- list(
  Category = c(
    "exacerbated"        = "#54278f",
    "no mitigation"      = "#469de1",
    "weak mitigation"    = "#ffd27f",
    "partial mitigation" = "#ff9e4e",
    "full mitigation"    = "#ff594e",
    "reversed"           = "#08519c"
  )
)

RdBu <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(30)
prox_to_dist_shift_breaks <- seq(-0.7, 0.7, length.out = 31)

# heatmap
pheatmap(
  heatmap_data,
  clustering_distance_cols = 'manhattan',
  clustering_distance_rows = "binary",
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = c("DMAi (mean)", "DMAi + PCF11 c1", "DMAi + PCF11 c2"),
  color  = RdBu,
  breaks = prox_to_dist_shift_breaks,
  fontsize_row = 1.25,
  annotation_row    = annotations,
  annotation_colors = annotation_colors,
  angle_col = 45
)

pheatmap(
  heatmap_data_with_gfp,
  clustering_distance_cols = 'manhattan',
  clustering_distance_rows = "binary",
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = c("DMAi (mean)", "DMAi + PCF11 c1", "DMAi + PCF11 c2", "DMAi + GFP"),
  color  = RdBu,
  breaks = prox_to_dist_shift_breaks,
  fontsize_row = 1.25,
  annotation_row    = annotations,
  annotation_colors = annotation_colors,
  gaps_col = c(1, 3),
  angle_col = 45
)

# =============================================================================
# SECTION 6: GENE GROUPS -> FASTA (site-level, two rows per gene)
# =============================================================================

strong_mit_genes    <- TUTR_categorized_both_rows %>% dplyr::filter(motif_group == "strong_mitigation")
not_mit_genes       <- TUTR_categorized_both_rows %>% dplyr::filter(motif_group == "not_mitigated")
exac_genes          <- TUTR_categorized_both_rows %>% dplyr::filter(mitigation_category == "exacerbated")
all_mitigated_genes <- TUTR_categorized_both_rows %>%
  dplyr::filter(mitigation_category %in% c("full mitigation", "reversed", "partial mitigation"))
dmai_lengthened_genes <- merged_df %>%
  dplyr::filter(gene_name %in% dmai_lengthened)

gene_group_list <- list(
  strong_mitigation = strong_mit_genes,
  not_mitigated     = not_mit_genes,
  all_mitigated     = all_mitigated_genes
)

#make fasta
genome <- BSgenome.Hsapiens.UCSC.hg38

extract_fasta <- function(bed_df, genome, out_file) {
  chrs <- bed_df$chr
  if (!all(grepl("^chr", chrs))) chrs <- paste0("chr", chrs)
  
  gr <- GenomicRanges::GRanges(
    seqnames = chrs,
    ranges   = IRanges::IRanges(start = bed_df$region_start + 1,  # BED 0-based -> GRanges 1-based
                                end   = bed_df$region_end),
    strand   = bed_df$strand
  )
  names(gr) <- paste0(bed_df$gene_name, "|", bed_df$feature_id)
  
  seqs <- BSgenome::getSeq(genome, gr)   # auto revcomp on minus strand
  names(seqs) <- names(gr)
  Biostrings::writeXStringSet(seqs, filepath = out_file)
  seqs
}

fasta_out_dir <- "Data/Figure_4_S7/fasta/1000nt_window"
dir.create(fasta_out_dir, recursive = TRUE, showWarnings = FALSE)

for (nm in names(all_APA_sites_with_window_list)) {
  df <- all_APA_sites_with_window_list[[nm]]
  if (nrow(df) == 0) next
  extract_fasta(df, genome,
                out_file = file.path(fasta_out_dir, paste0("PCF11_", nm, "_sites.fa")))
}

# =============================================================================
# HEADROOM AND INTER-SITE DISTANCE
#            (strong mitigation vs not mitigated; partial excluded)
# =============================================================================

two_group_df <- TUTR_categorized_both_rows %>%
  dplyr::filter(motif_group %in% c("strong_mitigation", "not_mitigated")) %>%
  dplyr::mutate(
    motif_group = factor(motif_group,
                         levels = c("not_mitigated", "strong_mitigation"),
                         labels = c("Not mitigated", "Strong mitigation"))
  )

group_cols <- c("Not mitigated" = "#ff594e", "Strong mitigation" = "#469de1")

# ---- 8a. Baseline proximal usage in DMAi (headroom control) ----
prox_usage_df <- two_group_df %>% dplyr::filter(pA_type == "proximal pA")
n_usage       <- prox_usage_df %>% dplyr::count(motif_group)

prox_usage_df %>%
  dplyr::group_by(motif_group) %>%
  dplyr::summarise(median_prox_usage = median(dmai_prox_usage, na.rm = TRUE),
                   n = dplyr::n()) %>%
  print()

ggplot(prox_usage_df, aes(x = motif_group, y = dmai_prox_usage, fill = motif_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.4, colour = "black") +
  geom_text(data = n_usage, aes(x = motif_group, y = 1.03, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 4.5) +
  ggpubr::stat_compare_means(method = "wilcox.test", label.y = 1.12, size = 5) +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(x = NULL,
       y = "Proximal pA usage in DMAi\n(no dox)",
       title = "Is PCF11 mitigation explained by headroom?") +
  theme_classic(base_size = 17) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16))

# ---- 8b. Proximal-distal inter-site distance ----
site_distance <- two_group_df %>%
  dplyr::mutate(
    pa_coord = dplyr::if_else(as.character(strand) == "+",
                              as.numeric(end), as.numeric(start))
  ) %>%
  dplyr::select(gene_name, motif_group, pA_type, pa_coord) %>%
  tidyr::pivot_wider(names_from = pA_type, values_from = pa_coord) %>%
  dplyr::rename(prox_coord = `proximal pA`, dist_coord = `distal pA`) %>%
  dplyr::mutate(site_distance = abs(dist_coord - prox_coord)) %>%
  dplyr::filter(!is.na(site_distance))

n_dist <- site_distance %>% dplyr::count(motif_group)

site_distance %>%
  dplyr::group_by(motif_group) %>%
  dplyr::summarise(median_distance = median(site_distance, na.rm = TRUE),
                   n = dplyr::n()) %>%
  print()

y_max <- max(site_distance$site_distance, na.rm = TRUE)

ggplot(site_distance, aes(x = motif_group, y = site_distance, fill = motif_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_text(data = n_dist, aes(x = motif_group, y = y_max * 0.8, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 4.5) +
  scale_fill_manual(values = group_cols) +
  scale_x_discrete(
    labels = c(
      "Not mitigated" = "Not PCF11-oe mitigated",
      "Strong mitigation" = "Strong PCF11-oe mitigation"
    )
  ) +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL,
       y = "pPA to dPA distance (bp)")+
  theme_classic(base_size = 17) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        ))

wilcox.test(
  site_distance ~ motif_group,
  data = site_distance
)

# ---- 8c. Continuous model: does headroom or spacing predict mitigation? ----
model_df <- prox_usage_df %>%
  dplyr::select(gene_name, ratio_mean, dmai_prox_usage, shift_alone_mean) %>%
  dplyr::inner_join(site_distance %>% dplyr::select(gene_name, site_distance),
                    by = "gene_name") %>%
  dplyr::mutate(log_distance = log10(site_distance))

mitigation_model <- lm(ratio_mean ~ dmai_prox_usage + log_distance + shift_alone_mean,
                       data = model_df)
print(summary(mitigation_model))

# =============================================================================
# SECTION 9: RECONCILIATION CHECK
# Group sizes in the BED source must match the plot groups.
# =============================================================================

cat("\nBED source (genes):\n")
print(sapply(gene_group_list[c("strong_mitigation", "not_mitigated")],
             function(df) dplyr::n_distinct(df$gene_name)))

cat("\nPlot groups (genes):\n")
print(two_group_df %>% dplyr::filter(pA_type == "proximal pA") %>% dplyr::count(motif_group))

