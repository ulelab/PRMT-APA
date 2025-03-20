#load libraries
library(megadepth)
library(recount3)
library(data.table)
library(optparse)


#######################################################
# 1) Parse Multiple Samples
#######################################################

#manually input samples
project_value <- "SRP017305"            # Example project
sample_ids <- c("SRR1238549", "SRR1238550", "SRR1238551", 
                "SRR1238552")  # Vector of sample IDs
output_dir <- "SRA_coverage_results/SRP017305/siCFIM25_SRP017305_results_whole_gene_approach_test/"  # Output folder

#######################################################
# 2) Coverage Calculation Function (returns full vector)
#######################################################

calculate_median_proximal_use <- function(
    bigwig_url,
    proximal_bed,
    distal_bed,
    coverage_density_threshold = 0.066  # equivalent to ~≥10 reads in 150 nt
) {
  # 1) Read the proximal BED
  prox_df <- data.table::fread(proximal_bed, header = FALSE)
  data.table::setnames(prox_df, c("chr", "start", "end", "gene_name", "length", "strand"))
  
  # 2) Read the distal BED
  dist_df <- data.table::fread(distal_bed, header = FALSE)
  data.table::setnames(dist_df, c("chr", "start", "end", "gene_name", "length", "strand"))
  
  # Check for unique gene names
  if (any(duplicated(prox_df$gene_name))) {
    stop("Duplicate gene names found in proximal BED file.")
  }
  if (any(duplicated(dist_df$gene_name))) {
    stop("Duplicate gene names found in distal BED file.")
  }
  
  # 3) Coverage for proximal and distal
  bw_cov_prox_gr <- megadepth::get_coverage(
    bigwig_file = bigwig_url,
    annotation  = proximal_bed,
    op = "sum"
  )
  bw_cov_dist_gr <- megadepth::get_coverage(
    bigwig_file = bigwig_url,
    annotation  = distal_bed,
    op = "sum"
  )
  
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
  
  # Instead of just returning numeric, we keep gene_name + ratio:
  all_proximal_use <- filtered_df[, .(gene_name, proximal_use)]
  
  return(list(
    median_proximal_use = median_prox_use,
    sites_quantified    = sites_count,
    all_proximal_use    = all_proximal_use
  ))
}

#######################################################
# 3) Condition -> BED Mappings
#######################################################
conditions <- list(
  DMAi = c("bed/first_exon_to_pPA/DMAi_sites_first_exon_to_pPA_sorted.bed", 
           "bed/UTR_regions/DMAi_sites_UTRregions_sorted.bed"),
  non_DMAi = c("bed/first_exon_to_pPA/non_DMAi_sites_first_exon_to_pPA_sorted.bed", 
               "bed/UTR_regions/non_DMAi_sites_UTRregions_sorted.bed"),
  control = c("bed/first_exon_to_pPA/control_sites_first_exon_to_pPA_sorted.bed", 
              "bed/UTR_regions/control_sites_UTRregions_sorted.bed"),
  siCFIM25 = c("bed/first_exon_to_pPA/siCFIM25_sites_first_exon_to_pPA_sorted.bed", 
               "bed/UTR_regions/siCFIM25_sites_UTRregions_sorted.bed"))

#######################################################
# 4) Main Loop: Collect Results, Then Subset
#######################################################

# We'll store results for each sample in memory
all_samples_results <- list()

# Process each single-SRR sample ID from `sample_ids`
for (sid in sample_ids) {
  
  cat("\n=== Processing SRR:", sid, "===\n")
  
  # Locate the BigWig URL for this sample
  bigwig_url <- recount3::locate_url(
      project_value,
      "data_sources/sra",
      "bw",
      "human",
      sid,
      recount3_url  = "https://recount-opendata.s3.amazonaws.com/recount3/release")
  
  # Prepare containers for summary + full results
  summary_list <- list()
  full_list    <- list()
  
  # For each condition, compute coverage
  for (cond in names(conditions)) {
    prox_bed <- conditions[[cond]][1]
    dist_bed <- conditions[[cond]][2]
    
    # Calculate usage
    res <- calculate_median_proximal_use(
      bigwig_url,
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
  
  # Store in memory (we won't write yet)
  all_samples_results[[sid]] <- list(
    summary_df = summary_df,
    full_df    = full_df
  )
}

#######################################################
# 5) Intersection of gene_names Across All SRRs
#######################################################
# Only keep gene_names that pass the coverage threshold in *every* sample.

# 1) Gather the gene_name sets from each sample's full_df
gene_sets <- lapply(all_samples_results, function(x) {
  unique(x$full_df$gene_name)
})

# 2) Compute the intersection
common_genes <- Reduce(intersect, gene_sets)

cat("\nNumber of common genes passing threshold in all samples:",
    length(common_genes), "\n")

#######################################################
# 6) Subset + Write Final Output
#######################################################
for (sid in sample_ids) {
  
  cat("\n=== Finalizing results for SRR:", sid, "===\n")
  
  # Retrieve the stored summary/full data
  summary_df <- all_samples_results[[sid]]$summary_df
  full_df    <- all_samples_results[[sid]]$full_df
  
  # Subset 'full_df' to only keep common_genes
  full_df <- full_df[ full_df$gene_name %in% common_genes, ]
  
  # (Optional) Recompute summary stats on this intersection
  # to reflect the new site count and median usage
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
  
  # Write the updated summary
  summary_file <- file.path(output_dir, paste0(sid, "_PA_usage_quant_summary.csv"))
  data.table::fwrite(summary_df_final, file = summary_file, quote = FALSE, 
                     sep = ",", col.names = TRUE)
  cat("  Updated summary written to:", summary_file, "\n")
  
  # Write the updated full output
  # We'll keep columns in a user-friendly order
  full_df <- full_df[, c("sample", "condition", "gene_name", "proximal_use")]
  
  full_file <- file.path(output_dir, paste0(sid, "_PA_usage_full_quant_results.csv"))
  data.table::fwrite(full_df, file = full_file, quote = FALSE, 
                     sep = ",", col.names = TRUE)
  cat("  Updated full results written to:", full_file, "\n")
}

cat("\nAll samples completed. Intersection-based filtering done.\n")



