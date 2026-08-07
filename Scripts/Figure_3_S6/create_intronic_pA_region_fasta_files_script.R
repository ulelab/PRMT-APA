library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(data.table)

#read in files
sig_dmai   <- read_csv("Data/Figure_1_S1_S2_S3_S4/APA_CSV_files/timecourse/APA_classified/sig_CSV/sig_DMAi_72hrs_iAPA.csv", show_col_types = FALSE)
ctrl_dmai  <- read_csv("Data/Figure_1_S1_S2_S3_S4/APA_CSV_files/timecourse/APA_classified/ctrl_CSV/ctrl_DMAi_72hrs_iAPA.csv", show_col_types = FALSE)

# ── Helper functions
collapse_gene_effect <- function(df, source_name) {
  df %>%
    dplyr::filter(!is.na(gene_name),
                  !is.na(internal_to_UTR_shift),
                  internal_to_UTR_shift != 0) %>%
    dplyr::distinct(gene_name, internal_to_UTR_shift) %>%
    dplyr::mutate(
      source     = source_name,
      abs_effect = abs(internal_to_UTR_shift),
      direction  = dplyr::case_when(
        internal_to_UTR_shift > 0 ~ "decreased",
        internal_to_UTR_shift < 0 ~ "increased"
      )
    ) %>%
    dplyr::select(gene_name, source, internal_to_UTR_shift, abs_effect, direction)
}

resolve_overlap <- function(df, threshold = 0.75) {
  if (nrow(df) == 1) {
    return(df %>%
             dplyr::mutate(resolution        = "unique",
                           competitor_effect = NA_real_,
                           competitor_ratio  = NA_real_))
  }
  df <- df %>%
    dplyr::arrange(dplyr::desc(abs_effect)) %>%
    dplyr::distinct(source, .keep_all = TRUE)
  if (nrow(df) == 1) {
    return(df %>%
             dplyr::mutate(resolution        = "unique_after_dedup",
                           competitor_effect = NA_real_,
                           competitor_ratio  = NA_real_))
  }
  df          <- df %>% dplyr::arrange(dplyr::desc(abs_effect))
  top_effect  <- df$abs_effect[1]
  sec_effect  <- df$abs_effect[2]
  ratio       <- sec_effect / top_effect
  if (ratio < threshold) {
    return(df %>%
             dplyr::slice(1) %>%
             dplyr::mutate(resolution        = "kept_stronger",
                           competitor_effect = sec_effect,
                           competitor_ratio  = ratio))
  } else {
    return(tibble(
      source                = character(),
      internal_to_UTR_shift = numeric(),
      abs_effect            = numeric(),
      direction             = character(),
      candidate_group       = character(),
      resolution            = character(),
      competitor_effect     = numeric(),
      competitor_ratio      = numeric()
    ))
  }
}

calculate_PA_regions_1000 <- function(df) {
  df %>%
    dplyr::mutate(
      pa_coord     = if_else(strand == "+", end, start),
      region_start = pmax(pa_coord - 500, 0),
      region_end   = pa_coord + 500
    )
}

prepare_bed_df <- function(df) {
  df %>%
    dplyr::transmute(
      chr        = chr,
      start      = region_start,
      end        = region_end,
      gene_name  = gene_name,
      feature_id = feature_id,
      strand     = strand
    ) %>%
    dplyr::distinct()
}

write_bed <- function(df, file) {
  write.table(df, file = file, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ── Create gene-level candidate groups
sig_dmai_gene   <- collapse_gene_effect(sig_dmai,   "DMAi")

candidate_tbl <- dplyr::bind_rows(sig_dmai_gene, sig_cdk11i_gene) %>%
  dplyr::mutate(candidate_group = paste(source, direction, "iAPA", sep = "_"))

# ── Resolve overlaps
assigned_treatment_tbl <- candidate_tbl %>%
  dplyr::group_by(gene_name) %>%
  dplyr::group_modify(~ resolve_overlap(.x, threshold = 0.75)) %>%
  dplyr::ungroup()

overlap_removed_tbl <- candidate_tbl %>%
  dplyr::group_by(gene_name) %>%
  dplyr::filter(dplyr::n_distinct(source) > 1) %>%
  dplyr::summarise(
    max_effect = max(abs_effect),
    min_effect = min(abs_effect),
    ratio      = min_effect / max_effect,
    .groups    = "drop"
  ) %>%
  dplyr::filter(ratio >= 0.75)

# ── Final gene sets
DMAi_increased_iAPA   <- assigned_treatment_tbl %>% dplyr::filter(candidate_group == "DMAi_increased_iAPA")   %>% dplyr::pull(gene_name)

# ── DMAi-specific unchanged: genes in the DMAi control set, excluding DMAi-significant genes
DMAi_sig_genes      <- unique(c(DMAi_increased_iAPA, DMAi_decreased_iAPA))
DMAi_unchanged_iAPA <- setdiff(unique(ctrl_dmai$gene_name), DMAi_sig_genes)

# ── Row-level data frames
DMAi_increased_rows   <- sig_dmai   %>% dplyr::filter(gene_name %in% DMAi_increased_iAPA)
DMAi_unchanged_rows <- ctrl_dmai %>% dplyr::filter(gene_name %in% DMAi_unchanged_iAPA)

groups_list <- list(
  "DMAi increased"   = DMAi_increased_rows,
  "DMAi unchanged" = DMAi_unchanged_rows
)

# ── FASTA extraction
genome <- BSgenome.Hsapiens.UCSC.hg38

extract_bed_to_fasta <- function(bed_df, genome, out_fasta) {
  gr <- GRanges(
    seqnames = bed_df$chr,
    ranges   = IRanges(start = bed_df$start + 1, end = bed_df$end),
    strand   = bed_df$strand
  )
  seqs <- getSeq(genome, gr)
  names(seqs) <- paste0(bed_df$feature_id, "|", bed_df$gene_name, "|",
                        bed_df$chr, ":", bed_df$start, "-", bed_df$end,
                        "|", bed_df$strand)
  writeXStringSet(seqs, filepath = out_fasta)
  seqs
}

setwd("Data/Figure_3_S6/fasta/1000nt_window/LU99_timecourse")

bed_list <- list(
  DMAi_increased_intronic   = DMAi_increased_rows   %>% dplyr::filter(pA_type == "proximal pA") %>% calculate_PA_regions_1000() %>% prepare_bed_df(),
  DMAi_increased_UTR        = DMAi_increased_rows   %>% dplyr::filter(pA_type == "distal pA")   %>% calculate_PA_regions_1000() %>% prepare_bed_df(),
  DMAi_unchanged_intronic = DMAi_unchanged_rows %>% dplyr::filter(pA_type == "proximal pA") %>% calculate_PA_regions_1000() %>% prepare_bed_df(),
  DMAi_unchanged_UTR      = DMAi_unchanged_rows %>% dplyr::filter(pA_type == "distal pA") %>% calculate_PA_regions_1000() %>% prepare_bed_df()
)

fasta_list <- lapply(names(bed_list), function(nm) {
  extract_bed_to_fasta(bed_list[[nm]], genome,
                       file.path("fasta", paste0(nm, ".fa")))
})
