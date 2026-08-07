# Load required packages
library(tidyverse)       # For data manipulation and visualization
library(tximport)        # For importing transcript-level estimates
library(DESeq2)          # For differential expression analysis
library(apeglm)          # For log fold change shrinkage
library(AnnotationDbi)   # For annotation functions
library(org.Hs.eg.db)    # For human gene annotations
library(GenomicFeatures) # For creating TxDb objects
library(ggplot2)         # For figures
library(RColorBrewer)    # For figures
library(pheatmap)        # For figures

# Read in sample metadata e.g
LU99_PCF11_sample_info <- read_csv('Data/Figure_4_S7/LU99_PCF11_DESeq2_sample_info.csv') %>%
  mutate(Treatment = factor(Treatment))

# Create vector of Salmon quant.sf file paths and assign sample names
LU99_PCF11_sample_files <- paste0(LU99_PCF11_sample_info$Treatment, LU99_PCF11_sample_info$Replicate, '/quant.sf')
names(LU99_PCF11_sample_files) <- paste0(LU99_PCF11_sample_info$Treatment, '_', LU99_PCF11_sample_info$Replicate)

# Base path relative to the *current working directory*
base_counts_dir <- "Data/Figure_4_S7/salmon_files_for_DESeq2/"

# Create vector of Salmon quant.sf file paths and assign sample names
LU99_PCF11_sample_files <- file.path(
  base_counts_dir,
  paste0(LU99_PCF11_sample_info$Treatment, LU99_PCF11_sample_info$Replicate),
  "quant.sf"
)
names(LU99_PCF11_sample_files) <- paste0(LU99_PCF11_sample_info$Treatment, '_', LU99_PCF11_sample_info$Replicate)

# Create tx2gene mapping using GTF annotation
txdb = txdbmaker::makeTxDbFromGFF(file='/path/to/gencode.v45.annotation.gtf')
k = keys(txdb, keytype = "GENEID")
df = select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene_map = df[, 2:1]
tx2gene_map$TXNAME = sub("\\..*", "", tx2gene_map$TXNAME)

# Import transcript-level estimates and summarize to gene-level counts
txi <- tximport(files= LU99_PCF11_sample_files,
                type = 'salmon',
                tx2gene = tx2gene_map,
                countsFromAbundance="no",
                ignoreTx = TRUE)

# Round counts to integers (DESeq2 requires integer counts)
txi$counts <- round(txi$counts)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = txi$counts,
                              colData = LU99_PCF11_sample_info,
                              design = ~ Treatment)

# Filter out low-count genes (keep genes with counts >5 in at least two samples)
keep <- rowSums(counts(dds) > 5) >= 2
dds <- dds[keep, ]

# Run DESeq normalization and estimation
dds <- DESeq(dds)

legend_order <- c(
  # No dox + DMSO
  "LU99_GFP_no_dox_DMSO",
  "LU99_PCF11_c3_no_dox_DMSO",
  "LU99_PCF11_c7_no_dox_DMSO",
  
  # + Dox + DMSO
  "LU99_GFP_5ng_dox_DMSO",
  "LU99_PCF11_c3_30ng_dox_DMSO",
  "LU99_PCF11_c7_20ng_dox_DMSO",
  
  # No dox + DMAi
  "LU99_GFP_no_dox_SA",
  "LU99_PCF11_c3_no_dox_SA",
  "LU99_PCF11_c7_no_dox_SA",
  
  # + Dox + DMAi
  "LU99_GFP_5ng_dox_SA",
  "LU99_PCF11_c3_30ng_dox_SA",
  "LU99_PCF11_c7_20ng_dox_SA"
)

#check PCF11 expression
label_condition <- function(x) {
  dplyr::case_when(
    grepl("no_dox_DMSO", x) ~ "No dox + DMSO",
    grepl("no_dox_SA",   x) ~ "No dox + DMAi",
    grepl("dox_DMSO",    x) ~ "Dox + DMSO",     # no_dox_* already matched above
    grepl("dox_SA",      x) ~ "Dox + DMAi",
    TRUE                    ~ NA_character_
  )
}

condition_levels <- c("No dox + DMSO", "No dox + DMAi",
                      "Dox + DMSO",    "Dox + DMAi")

condition_cols <- c(
  "No dox + DMSO" = "#BDBDBD",
  "No dox + DMAi" = "#CB181D",
  "Dox + DMSO"    = "#3B82F6",
  "Dox + DMAi"    = "#7E22CE"
)

# =============================================================================
# 2. Transgene counts (bowtie2, construct-specific references)
# =============================================================================

transgene_raw <- dplyr::bind_rows(
  readr::read_tsv("Data/Figure_4_S7/flag_pcf11_raw_counts.tsv",
                  col_names = c("sample_id", "count"), show_col_types = FALSE),
  readr::read_tsv("Data/Figure_4_S7/gfp_raw_counts.tsv",
                  col_names = c("sample_id", "count"), show_col_types = FALSE)
)

sf_df <- tibble::tibble(
  sample_id   = names(DESeq2::sizeFactors(dds)),
  size_factor = as.numeric(DESeq2::sizeFactors(dds))
) %>%
  dplyr::mutate(sample_id = gsub("__", "_", sample_id))

meta_df <- LU99_PCF11_sample_info %>%
  dplyr::mutate(sample_id = gsub("__", "_", paste0(Treatment, "_", Replicate))) %>%
  dplyr::select(sample_id, Treatment)

transgene <- transgene_raw %>%
  dplyr::left_join(sf_df,   by = "sample_id") %>%
  dplyr::left_join(meta_df, by = "sample_id") %>%
  dplyr::mutate(
    norm_count = count / size_factor,
    Condition  = label_condition(Treatment),
    Line = dplyr::case_when(
      grepl("GFP",      Treatment) ~ "GFP",
      grepl("PCF11_c3", Treatment) ~ "PCF11 c1",
      grepl("PCF11_c7", Treatment) ~ "PCF11 c2"
    ),
    Condition = factor(Condition, levels = condition_levels),
    Line      = factor(Line, levels = c("GFP", "PCF11 c1", "PCF11 c2"))
  )

# ---- Induction / leak check ----
print(
  transgene %>%
    dplyr::group_by(Line, Condition) %>%
    dplyr::summarise(median = median(norm_count), min = min(norm_count),
                     max = max(norm_count), n = dplyr::n(), .groups = "drop")
)

# =============================================================================
# 3. Transgene expression figure
# =============================================================================

p_transgene <- ggplot(transgene,
                      aes(x = Condition, y = norm_count, colour = Condition)) +
  geom_point(size = 2.5, position = position_jitter(width = 0.12, height = 0)) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, linewidth = 0.4, colour = "black") +
  facet_wrap(~ Line, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = condition_cols) +
  labs(x = NULL,
       y = "Transgene counts\n(size factor-normalised)") +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 13),
    strip.text   = element_text(size = 17, face = "bold"),
    plot.caption = element_text(size = 9, hjust = 0),
    legend.position = "none"
  )
print(p_transgene)

# ---- Fold change: DMAi vs DMSO within dox-induced samples ----
fold_changes <- transgene %>%
  dplyr::filter(Condition %in% c("Dox + DMSO", "Dox + DMAi")) %>%
  dplyr::group_by(Line, Condition) %>%
  dplyr::summarise(mean_norm = mean(norm_count), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Condition, values_from = mean_norm) %>%
  dplyr::mutate(fold_change = `Dox + DMAi` / `Dox + DMSO`)
print(fold_changes)

# =============================================================================
# Fold change: transgene expression, DMAi vs DMSO within dox-induced samples
# =============================================================================

# ---- Per-replicate ratios ----
# Replicates aren't paired across drug conditions, so pair them by replicate
# number. If replicate identity is not meaningful (independent wells), the
# alternative below (ratio of means with propagated error) is more honest.

dox_only <- transgene %>%
  dplyr::filter(Condition %in% c("Dox + DMSO", "Dox + DMAi")) %>%
  dplyr::mutate(
    replicate = stringr::str_extract(sample_id, "\\d+$")
  )

# ---- Option A: ratio of group means, with bootstrap CI ----
fc_summary <- dox_only %>%
  dplyr::group_by(Line, Condition) %>%
  dplyr::summarise(mean_norm = mean(norm_count),
                   sd_norm   = sd(norm_count),
                   n         = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Condition,
                     values_from = c(mean_norm, sd_norm, n)) %>%
  dplyr::mutate(
    fold_change = `mean_norm_Dox + DMAi` / `mean_norm_Dox + DMSO`,
    # error propagation for a ratio: relative errors add in quadrature
    rel_se_dmai = (`sd_norm_Dox + DMAi`  / sqrt(`n_Dox + DMAi`))  / `mean_norm_Dox + DMAi`,
    rel_se_dmso = (`sd_norm_Dox + DMSO`  / sqrt(`n_Dox + DMSO`))  / `mean_norm_Dox + DMSO`,
    se_fc       = fold_change * sqrt(rel_se_dmai^2 + rel_se_dmso^2),
    ymin        = fold_change - se_fc,
    ymax        = fold_change + se_fc
  )

print(fc_summary %>% dplyr::select(Line, fold_change, se_fc))

# ---- Individual replicate ratios, for overlaying as points ----
# Each DMAi replicate divided by the MEAN of the DMSO replicates for that line
dmso_means <- dox_only %>%
  dplyr::filter(Condition == "Dox + DMSO") %>%
  dplyr::group_by(Line) %>%
  dplyr::summarise(dmso_mean = mean(norm_count), .groups = "drop")

fc_points <- dox_only %>%
  dplyr::filter(Condition == "Dox + DMAi") %>%
  dplyr::left_join(dmso_means, by = "Line") %>%
  dplyr::mutate(fold_change = norm_count / dmso_mean)

# =============================================================================
# Figure
# =============================================================================

line_cols <- c("GFP" = "#35BD38", "PCF11 c1" = "#7E22CE", "PCF11 c2" = "#9D4EDD")

p_fc <- ggplot(fc_summary, aes(x = Line, y = fold_change, fill = Line)) +
  geom_col(width = 0.65, alpha = 0.85, colour = "black", linewidth = 1) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, linewidth = 0.8) +
  geom_point(data = fc_points, aes(x = Line, y = fold_change),
             inherit.aes = FALSE, size = 2.5, colour = "black",
             position = position_jitter(width = 0.08, height = 0)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40")  +
  scale_fill_manual(values = line_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL,
       y = "Transgene expression\n(DMAi / DMSO, dox-induced)") +
  theme_classic(base_size = 17) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 15),
    plot.caption = element_text(size = 9, hjust = 0)
  )
print(p_fc)