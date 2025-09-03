#Author: Charlotte Capitanchik

library(megadepth)
library(recount3)
library(data.table)
library(optparse)

#setwd
setwd("/Users/k2362866/Documents/AZ_postdoc/files_for_charlotte/github/motif_and_clip_analysis/CFIM_mitigation_categories/")

# Bigwig file path
bigwig_base <- "/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/Charlotte_metaprofiles/RNA_seq_files/ctrl"

sample_bigwig_map <- list(
  "shCtrl-K562-R1" = c(
    file.path(bigwig_base, "shCtrl_K562_1_unique_pos.bigWig"),
    file.path(bigwig_base, "shCtrl_K562_1_unique_neg.bigWig")
  ),
  "shCtrl-K562-R2" = c(
    file.path(bigwig_base, "shCtrl_K562_2_unique_pos.bigWig"),
    file.path(bigwig_base, "shCtrl_K562_2_unique_neg.bigWig")
  ),
  "shCtrl-HepG2-R1" = c(
    file.path(bigwig_base, "shCtrl_HepG2_1_unique_pos.bigWig"),
    file.path(bigwig_base, "shCtrl_HepG2_1_unique_neg.bigWig")
  ),
  "shCtrl-HepG2-R2" = c(
    file.path(bigwig_base, "shCtrl_HepG2_2_unique_pos.bigWig"),
    file.path(bigwig_base, "shCtrl_HepG2_2_unique_neg.bigWig")
  )
)

#######################################################
# 2) Coverage Calculation Function
#######################################################

get_combined_coverage <- function(bw_files, annotation_bed, op = "sum") {
  # For each BigWig in 'bw_files', compute coverage with 'megadepth::get_coverage()'
  coverage_list <- lapply(bw_files, function(bw) {
    megadepth::get_coverage(
      bigwig_file = bw,
      annotation  = annotation_bed,
      op          = op
    )
  })
  
  # Start with coverage from the first BigWig
  combined <- coverage_list[[1]]
  
  # Sum coverage from any additional BigWigs
  if (length(coverage_list) > 1) {
    for (i in seq(2, length(coverage_list))) {
      combined$score <- combined$score + coverage_list[[i]]$score
    }
  }
  
  return(combined)
}


#######################################################
# 3) Condition -> BED Mappings
#######################################################
conditions <- list(
  DMAi_full_mit = c("bed/first_exon_to_pPA/expression_matched/DMAi_full_mitigation_sites_first_exon_to_pPA_sorted.bed", 
                    "bed/UTR_regions/expression_matched/DMAi_full_mitigation_sites_UTRregions_sorted.bed"),
  DMAi_part_mit = c("bed/first_exon_to_pPA/expression_matched/DMAi_part_mitigation_sites_first_exon_to_pPA_sorted.bed", 
                    "bed/UTR_regions/expression_matched/DMAi_part_mitigation_sites_UTRregions_sorted.bed"),
  DMAi_no_mit = c("bed/first_exon_to_pPA/expression_matched/DMAi_no_mitigation_sites_first_exon_to_pPA_sorted.bed", 
                  "bed/UTR_regions/expression_matched/DMAi_no_mitigation_sites_UTRregions_sorted.bed"))

#######################################################
# 4) Main Loop: We'll Collect Results, Then Subset
#######################################################

calculate_median_proximal_use <- function(
    bigwig_files,       # now a vector of local .bigwig paths
    proximal_bed,
    distal_bed,
    coverage_density_threshold = 0.066
) {
  # 1) Read the proximal BED
  prox_df <- data.table::fread(proximal_bed, header = FALSE)
  data.table::setnames(prox_df, c("chr", "start", "end", "gene_name", "length", "strand"))
  
  # 2) Read the distal BED
  dist_df <- data.table::fread(distal_bed, header = FALSE)
  data.table::setnames(dist_df, c("chr", "start", "end", "gene_name", "length", "strand"))
  
  # Check for duplicates
  if (any(duplicated(prox_df$gene_name))) {
    stop("Duplicate gene names found in proximal BED file.")
  }
  if (any(duplicated(dist_df$gene_name))) {
    stop("Duplicate gene names found in distal BED file.")
  }
  
  # 3) Coverage for proximal and distal, combining multiple BigWigs
  bw_cov_prox_gr <- get_combined_coverage(bigwig_files, proximal_bed, op = "sum")
  bw_cov_dist_gr <- get_combined_coverage(bigwig_files, distal_bed,   op = "sum")
  
  # 4) Convert GRanges -> data.table
  df_bw_cov_prox <- data.table::as.data.table(as.data.frame(bw_cov_prox_gr))
  data.table::setnames(
    df_bw_cov_prox,
    old = c("seqnames", "start", "end", "strand", "score"),
    new = c("chr", "start", "end", "strand", "coverage_prox")
  )
  
  df_bw_cov_dist <- data.table::as.data.table(as.data.frame(bw_cov_dist_gr))
  data.table::setnames(
    df_bw_cov_dist,
    old = c("seqnames", "start", "end", "strand", "score"),
    new = c("chr", "start", "end", "strand", "coverage_dist")
  )
  
  # 5) Attach coverage data to original BED data
  prox_df[, coverage_prox := df_bw_cov_prox$coverage_prox]
  dist_df[, coverage_dist := df_bw_cov_dist$coverage_dist]
  
  # 6) Compute coverage density
  prox_df[, coverage_density_prox := coverage_prox / length]
  dist_df[, coverage_density_dist := coverage_dist / length]
  
  # 7) Merge by gene_name
  merged_df <- merge(
    prox_df[, .(gene_name, coverage_density_prox)],
    dist_df[, .(gene_name, coverage_density_dist)],
    by = "gene_name",
    all = FALSE
  )
  
  # 8) Filter intervals by coverage density threshold
  filtered_df <- merged_df[
    coverage_density_prox >= coverage_density_threshold |
      coverage_density_dist >= coverage_density_threshold
  ]
  
  # 9) Compute ratio (avoid zero coverage by substituting small value)
  filtered_df[, proximal_use := fifelse(
    coverage_density_prox == 0, 0.001, coverage_density_prox
  )]
  filtered_df[, distal_use := fifelse(
    coverage_density_dist == 0, 0.001, coverage_density_dist
  )]
  
  filtered_df[, proximal_use := proximal_use / (proximal_use + distal_use)]
  
  # 10) Summaries
  median_prox_use <- median(filtered_df$proximal_use, na.rm = TRUE)
  sites_count     <- nrow(filtered_df)
  
  # Return gene_name + ratio
  all_proximal_use <- filtered_df[, .(gene_name, proximal_use)]
  
  return(list(
    median_proximal_use = median_prox_use,
    sites_quantified    = sites_count,
    all_proximal_use    = all_proximal_use
  ))
}

# 1) List your samples
sample_ids <- c("shCtrl-K562-R1","shCtrl-K562-R2","shCtrl-HepG2-R1","shCtrl-HepG2-R2")

# 2) We already have sample_bigwig_map for them
output_dir <- "SRA_coverage_results/shCtrl_K562_HepG2_CFIM_mitigation_categories"

# We'll store results for each sample in memory
all_samples_results <- list()

# 3) Loop over each sample
for (sid in sample_ids) {
  
  cat("\n=== Processing sample:", sid, "===\n")
  
  # Gather the local BigWig file paths
  bw_files <- sample_bigwig_map[[sid]]
  
  # Prepare containers for summary + full results
  summary_list <- list()
  full_list    <- list()
  
  # For each condition, compute coverage
  for (cond in names(conditions)) {
    prox_bed <- conditions[[cond]][1]
    dist_bed <- conditions[[cond]][2]
    
    # Calculate usage, passing 'bw_files' now
    res <- calculate_median_proximal_use(
      bigwig_files = bw_files,
      proximal_bed = prox_bed,
      distal_bed   = dist_bed
    )
    
    # Build summary
    summary_list[[cond]] <- data.frame(
      condition           = cond,
      median_proximal_use = round(res$median_proximal_use, 3),
      sites_quantified    = res$sites_quantified
    )
    
    # Build full data (includes gene_name)
    df_full <- data.frame(
      condition    = cond,
      gene_name    = res$all_proximal_use$gene_name,
      proximal_use = res$all_proximal_use$proximal_use
    )
    full_list[[cond]] <- df_full
  }
  
  # Combine per-sample results
  summary_df <- do.call(rbind, summary_list)
  summary_df$sample <- sid
  
  full_df <- do.call(rbind, full_list)
  full_df$sample <- sid
  
  # Store in memory
  all_samples_results[[sid]] <- list(
    summary_df = summary_df,
    full_df    = full_df
  )
}

# 4) Intersection of gene_names across all samples
gene_sets <- lapply(all_samples_results, function(x) {
  unique(x$full_df$gene_name)
})
common_genes <- Reduce(intersect, gene_sets)
cat("\nNumber of common genes passing threshold in all samples:", length(common_genes), "\n")

# 5) Subset + Write final output
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (sid in sample_ids) {
  cat("\n=== Finalizing results for sample:", sid, "===\n")
  
  summary_df <- all_samples_results[[sid]]$summary_df
  full_df    <- all_samples_results[[sid]]$full_df
  
  # Subset 'full_df' to only keep genes common across all samples
  full_df <- full_df[full_df$gene_name %in% common_genes, ]
  
  # Recompute summary for these common genes
  new_summary_list <- list()
  for (cond in unique(full_df$condition)) {
    sub_df  <- full_df[full_df$condition == cond, ]
    med_val <- median(sub_df$proximal_use, na.rm=TRUE)
    new_summary_list[[cond]] <- data.frame(
      sample             = sid,
      condition          = cond,
      median_proximal_use= round(med_val, 3),
      sites_quantified   = nrow(sub_df)
    )
  }
  summary_df_final <- do.call(rbind, new_summary_list)
  
  # Write updated summary
  summary_file <- file.path(output_dir, paste0(sid, "_PA_usage_quant_summary.csv"))
  data.table::fwrite(summary_df_final, file = summary_file,
                     quote = FALSE, sep = ",", col.names = TRUE)
  
  # Write updated full
  full_df <- full_df[, c("sample", "condition", "gene_name", "proximal_use")]
  full_file <- file.path(output_dir, paste0(sid, "_PA_usage_full_quant_results.csv"))
  data.table::fwrite(full_df, file = full_file,
                     quote = FALSE, sep = ",", col.names = TRUE)
}
cat("\nAll samples completed. Intersection-based filtering done.\n")
