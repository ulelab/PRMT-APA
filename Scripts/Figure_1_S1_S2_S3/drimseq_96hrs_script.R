# Load necessary libraries
library(data.table)     # For fast data manipulation
library(DRIMSeq)        # For differential transcript usage analysis
library(BiocParallel)   # For parallel computation
library(ggplot2)        # For plotting
library(ggpubr)         # For enhanced ggplot2 plots
library(tidyverse)      # For data manipulation and visualization
library(RColorBrewer)   # For color palettes
library(stageR)         # For stage-wise analysis in high-dimensional data

# Set working directory (change this to your own directory)
# Make sure all input files are located in this directory or provide full paths
setwd("/path/to/your/working/directory")
# Example:
setwd("/ulelab/PRMT-APA/Scripts/Figure_1_S1_S2_S3/")

# Step 1: Read in the PolyA Site (PAS) atlas file
# Input: PAS atlas file in BED format (e.g., 'merged_polya.filteredunique.annotated.bed')
# Ensure that this file contains the necessary columns as specified below
pas_atlas_file <- '../../Data/Figure_1_S1_S2_S3/pA_atlas/merged_polya.filteredunique.annotated.bed'
pas_data <- fread(pas_atlas_file, 
                       col.names = c("seqnames", "start", "end", "id", "score", "strand", "ensg", "hgnc", "region"))

# Calculate total score per gene (ensg)
pas_data[, total_score := sum(score), by = ensg]

# Calculate thresholds (5% and 1% of total score)
pas_data[, percent_5 := 0.05 * total_score]
pas_data[, percent_1 := 0.01 * total_score]

# Count number of PAS per gene
pas_data[, num_pas := .N, by = ensg]

# Step 2: Load and merge all count files
# Inputs:
# - Directory containing count files (e.g., 'counts/')
# - Pattern to match count files (e.g., '.bed')
counts_dir <- "../../Data/Figure_1_S1_S2_S3/counts/96hrs/"  # Directory containing count files
counts_pattern <- ".bed"  # Pattern to match count files (adjust if necessary)

# List all count files in the counts directory
all_count_files <- list.files(path = counts_dir, full.names = TRUE, pattern = counts_pattern)

# Read and merge count data from all samples
count_list <- lapply(all_count_files, function(file) {
  # Extract sample name from the filename
  # Adjust the pattern in gsub as needed to match your filenames
  sample_name <- gsub("_r1.polya_trimmed_window_merged_polya.bed", "", basename(file))
  # Read in the count data for this sample
  # Adjust 'select' indices and 'col.names' as per your file structure
  data <- fread(file, select = c(4, 7, 8, 10), 
                col.names = c("PAS", "ensg", "hgnc", sample_name))
  return(data)
})

# Merge all count data into a single data table by 'PAS', 'ensg', 'hgnc'
count_data <- Reduce(function(x, y) merge(x, y, by = c("PAS", "ensg", "hgnc")), count_list)

# Remove intergenic PAS
count_data <- count_data[ensg != "intergenic"]

# Keep only PAS that are present in genes with more than one PAS
count_data <- count_data[PAS %in% pas_data[num_pas > 1]$id]

# Create 'gene_id' and 'feature_id' columns required by DRIMSeq
count_data[, gene_id := paste0(hgnc, "_", ensg)]
count_data[, feature_id := PAS]

# Convert count data to data frame and remove unnecessary columns
# Remove 'PAS', 'ensg', 'hgnc' columns
counts_df <- as.data.frame(count_data[, -c(1:3)])

# Step 3: Read in sample metadata
# Input: Metadata file containing sample IDs and their corresponding conditions/groups
# Ensure the metadata file has columns 'SampleName' and 'Condition'
metadata_file <- "../../Data/Figure_1_S1_S2_S3/tso_seq_metadata.txt"
metadata <- read.table(metadata_file, header = TRUE, stringsAsFactors = FALSE)

# Create a data frame with sample IDs and their corresponding groups (conditions)
samples <- data.frame(sample_id = metadata$SampleName,
                      group = metadata$Condition)

# Step 4: Define your conditions and their corresponding groups
# Adjust the 'conditions_list' as per your dataset
# Each element in the list represents a condition with its control and treatment groups
conditions_list <- list(
  'ADMAi' = c('DMSO', 'ADMAi'),
  'SDMAi' = c('DMSO', 'SDMAi'),
  'SDMAi_ADMAi' = c('DMSO', 'SDMAi_ADMAi')
)

# Step 5: Create a function to process each condition
# This function encapsulates the per-condition processing steps
process_condition <- function(condition_name, groups, counts_df, samples) {
  # Subset samples for the condition
  samples_cond <- samples[samples$group %in% groups, ]
  samples_cond$group <- factor(samples_cond$group, levels = groups)
  
  # Subset counts data for the samples
  counts_df_cond <- counts_df[, c(samples_cond$sample_id, "gene_id", "feature_id")]
  
  # Create DRIMSeq data object
  d_cond <- dmDSdata(counts = counts_df_cond, samples = samples_cond)
  
  # Filter lowly expressed genes and features
  # Adjust the thresholds as per your requirements
  min_samps_gene_expr <- floor(nrow(samples_cond) * 0.75)
  min_samps_feature_expr <- floor(min(table(samples_cond$group)) * 0.75)
  d_cond <- dmFilter(d_cond, 
                     min_samps_gene_expr = min_samps_gene_expr, 
                     min_samps_feature_expr = min_samps_feature_expr, 
                     min_gene_expr = 10, 
                     min_feature_expr = 5)
  
  # Fit Dirichlet-multinomial models
  design_full_cond <- model.matrix(~ group, data = DRIMSeq::samples(d_cond))
  set.seed(42)  # For reproducibility
  d_cond <- dmPrecision(d_cond, design = design_full_cond, verbose = 1)
  d_cond <- dmFit(d_cond, design = design_full_cond, verbose = 1)
  
  # Test for differential transcript usage
  coef_cond <- colnames(design_full_cond)[2]
  d_cond <- dmTest(d_cond, coef = coef_cond, verbose = 1)
  
  # Handle NA p-values
  res_genes <- results(d_cond)
  res_pA_sites <- results(d_cond, level = "feature")
  res_genes$pvalue <- ifelse(is.na(res_genes$pvalue), 1, res_genes$pvalue)
  res_pA_sites$pvalue <- ifelse(is.na(res_pA_sites$pvalue), 1, res_pA_sites$pvalue)
  
  # Perform two-stage test using stageR
  pScreen <- res_genes$pvalue
  names(pScreen) <- res_genes$gene_id
  
  pConfirmation <- matrix(res_pA_sites$pvalue, ncol = 1)
  rownames(pConfirmation) <- res_pA_sites$feature_id
  
  tx2gene <- res_pA_sites[, c("feature_id", "gene_id")]
  
  stageR_obj <- stageRTx(pScreen = pScreen, 
                         pConfirmation = pConfirmation, 
                         pScreenAdjusted = FALSE, 
                         tx2gene = tx2gene)
  
  stageR_obj <- stageWiseAdjustment(object = stageR_obj, method = "dtu", alpha = 0.05)
  
  two_stage_padj <- getAdjustedPValues(stageR_obj, order = FALSE, onlySignificantGenes = FALSE)
  
  colnames(two_stage_padj) <- c('gene_id', 'feature_id', 'twostep_gene_padj', 'twostep_feature_padj')
  
  # Calculate change in feature usage
  proportions_cond <- proportions(d_cond)
  
  # Define control and treatment sample IDs
  control_samples <- samples_cond$sample_id[samples_cond$group == groups[1]]
  treatment_samples <- samples_cond$sample_id[samples_cond$group == groups[2]]
  
  # Calculate mean usage per group
  proportions_cond$pA_usage_control <- rowMeans(proportions_cond[, control_samples], na.rm = TRUE)
  proportions_cond$pA_usage_treatment <- rowMeans(proportions_cond[, treatment_samples], na.rm = TRUE)
  
  # Calculate change in usage
  proportions_cond$change_in_usage <- proportions_cond$pA_usage_treatment - proportions_cond$pA_usage_control
  
  # Merge results
  res_pA_sites_with_change <- left_join(
    res_pA_sites,
    proportions_cond[, c("gene_id", "feature_id", "pA_usage_control", "pA_usage_treatment", "change_in_usage")],
    by = c("gene_id", "feature_id")
  )
  
  # Combine two-stage adjusted p-values with change in usage
  two_stage_results <- inner_join(two_stage_padj, res_pA_sites_with_change, by = c('gene_id', 'feature_id'))
  
  # Identify significant features (adjusted p-value <= 0.05 and absolute change in usage >= 0.1)
  sig_two_stage <- two_stage_results %>% filter(twostep_feature_padj <= 0.05 & abs(change_in_usage) >= 0.1)
  
  # Return results as a list
  return(list(
    res_genes = res_genes,
    res_pA_sites = res_pA_sites,
    two_stage_padj = two_stage_padj,
    two_stage_results = two_stage_results,
    sig_two_stage = sig_two_stage,
    proportions_cond = proportions_cond,
    samples_cond = samples_cond,
    d_cond = d_cond,
    condition_name = condition_name
  ))
}

# Step 6: Process all conditions using the function
condition_results <- lapply(names(conditions_list), function(cond_name) {
  groups <- conditions_list[[cond_name]]
  result <- process_condition(cond_name, groups, counts_df, samples)
  return(result)
})
names(condition_results) <- names(conditions_list)

# Step 7: process atlas file
colnames(full_pas_data)[1:4] <- c("chr", "start", "end", "feature_id")
colnames(full_pas_data)[6] <- "strand"
polyA_coords <- full_pas_data[, c(1:4, 6)]

# Step 8: Merge results with PolyA coordinates and separate gene_id
for (i in seq_along(condition_results)) {
  cond_name <- condition_results[[i]]$condition_name
  
  # Merge two_stage_results with PolyA coordinates
  condition_results[[i]]$two_stage_results_merged <- inner_join(condition_results[[i]]$two_stage_results, polyA_coords, by = "feature_id")
  
  # Merge sig_two_stage with PolyA coordinates
  condition_results[[i]]$sig_two_stage_merged <- inner_join(condition_results[[i]]$sig_two_stage, polyA_coords, by = "feature_id")
  
  # Separate 'gene_id' into 'gene_name' and 'ensembl_ID'
  condition_results[[i]]$two_stage_results_merged <- condition_results[[i]]$two_stage_results_merged %>%
    separate(gene_id, into = c('gene_name', 'ensembl_ID'), sep = "_")
  condition_results[[i]]$sig_two_stage_merged <- condition_results[[i]]$sig_two_stage_merged %>%
    separate(gene_id, into = c('gene_name', 'ensembl_ID'), sep = "_")
}

# Step 9: Write results to files
# Output: CSV files containing all polyA sites and significant polyA sites with change in usage
# The files are saved in 'CSV_files/' directory
# Create output directories if they don't exist
dir.create("../../Data/Figure_1_S1_S2_S3/APA_CSV_files/96hrs/", showWarnings = FALSE)
dir.create("../../Data/Figure_1_S1_S2_S3/APA_CSV_files/96hrs/positional_info", showWarnings = FALSE)

for (i in seq_along(condition_results)) {
  cond_name <- condition_results[[i]]$condition_name
  
  # Define output file paths
  all_sites_file <- paste0("../../Data/Figure_1_S1_S2_S3/APA_CSV_files/96hrs/all_polyA_sites_", cond_name, "_with_change_in_usage.csv")
  sig_sites_file <- paste0("../../Data/Figure_1_S1_S2_S3/APA_CSV_files/96hrs/sig_polyA_sites_", cond_name, "_with_change_in_usage.csv")
  
  # Write to CSV
  write_csv(condition_results[[i]]$two_stage_results_merged, all_sites_file)
  write_csv(condition_results[[i]]$sig_two_stage_merged, sig_sites_file)
}

# Step 10: Identify top two most used PAS per gene and add positional info
# Function to identify top two PAS per gene based on mean usage
identify_top2_pas <- function(df) {
  # Calculate mean usage across control and treatment
  df$mean_usage <- (df$pA_usage_control + df$pA_usage_treatment) / 2
  # Group by gene and select top two PAS based on mean usage
  top2_pas <- df %>% group_by(gene_name) %>% slice_max(order_by = mean_usage, n=2)
  return(top2_pas)
}

# Function to identify positive strand transcripts
id_pos_top2 <- function(top2_df) {
  pos_transcripts <- top2_df %>% filter(strand == "+")
  return(pos_transcripts)
}

# Function to identify negative strand transcripts
id_neg_top2 <- function(top2_df) {
  neg_transcripts <- top2_df %>% filter(strand == "-")
  return(neg_transcripts)
}

# Function to find distal PAS in positive and negative strand transcripts and combine them
id_distal <- function(pos_df, neg_df) {
  pos_distal <- pos_df %>% group_by(gene_name) %>% slice_max(order_by = start)
  neg_distal <- neg_df %>% group_by(gene_name) %>% slice_min(order_by = start)
  distal_pAs <- bind_rows(pos_distal, neg_distal)
  return(distal_pAs)
}

# Function to find proximal PAS in positive and negative strand transcripts and combine them
id_proximal <- function(pos_df, neg_df) {
  pos_proximal <- pos_df %>% group_by(gene_name) %>% slice_min(order_by = start)
  neg_proximal <- neg_df %>% group_by(gene_name) %>% slice_max(order_by = start)
  proximal_pAs <- bind_rows(pos_proximal, neg_proximal)
  return(proximal_pAs)
}

# Function to add 'pA_type' column and combine distal and proximal PAS
unite_pAs <- function(dist_df, prox_df) {
  dist_df$pA_type <- "distal pA"
  prox_df$pA_type <- "proximal pA"
  all_pAs <- bind_rows(dist_df, prox_df)
  return(all_pAs)
}

# Process each condition to identify top two PAS and positional info
for (i in seq_along(condition_results)) {
  cond_name <- condition_results[[i]]$condition_name
  
  # Identify top two most used PAS per gene
  top2_pas <- identify_top2_pas(condition_results[[i]]$two_stage_results_merged)
  
  # Separate positive and negative strand transcripts
  pos_transcripts <- id_pos_top2(top2_pas)
  neg_transcripts <- id_neg_top2(top2_pas)
  
  # Identify distal and proximal PAS
  distal_pAs <- id_distal(pos_transcripts, neg_transcripts)
  proximal_pAs <- id_proximal(pos_transcripts, neg_transcripts)
  
  # Combine distal and proximal PAS and add 'pA_type' column
  combined_pAs <- unite_pAs(distal_pAs, proximal_pAs)
  
  # Add condition name for later use
  combined_pAs$condition <- cond_name
  
  # Store in the condition_results list
  condition_results[[i]]$combined_pAs <- combined_pAs
}

# Step 11: Calculate proximal to distal usage shift
# Function to calculate the shift between proximal and distal PAS usage
calculate_shift <- function(df) {
  df <- df %>%
    group_by(gene_name) %>%
    mutate(prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] - change_in_usage[pA_type == "proximal pA"]) %>%
    ungroup()
  return(df)
}

# Apply the shift calculation to each condition
for (i in seq_along(condition_results)) {
  condition_results[[i]]$shifted_pAs <- calculate_shift(condition_results[[i]]$combined_pAs)
}

# Step 12: Write positional info to files
# Output: CSV files containing top two PAS with positional information
for (i in seq_along(condition_results)) {
  cond_name <- condition_results[[i]]$condition_name
  output_file <- paste0("../../Data/Figure_1_S1_S2_S3/APA_CSV_files/96hrs/positional_info/top2_", cond_name, "_sites_with_positional_info.csv")
  write_csv(condition_results[[i]]$shifted_pAs, output_file)
}

# Step 13: Combine data for plotting
# Combine the shifted PAS data across conditions for plotting
combined_shifted_pas <- bind_rows(lapply(condition_results, function(res) res$shifted_pAs))

# Set factor levels for plotting
combined_shifted_pas$condition <- factor(combined_shifted_pas$condition, levels = c('ADMAi', 'SDMAi' ,'SDMAi_ADMAi'))
combined_shifted_pas$pA_type <- factor(combined_shifted_pas$pA_type, levels = c('proximal pA', 'distal pA'))

# Step 24: Plot the proximal to distal usage shift
# Remove top 5% highest absolute shifts in each group
filtered_data <- combined_shifted_pas %>%
  group_by(condition) %>%
  mutate(abs_shift = abs(prox_to_dist_shift)) %>%
  arrange(desc(abs_shift)) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(rank > 0.1 * max(rank)) %>%
  ungroup() %>%
  dplyr::select(-abs_shift, -rank) %>% # Remove temporary columns
  group_by(gene_name,condition) %>% dplyr::slice(1) %>% ungroup() # only need one row per gene

#define colour scheme
top2_colours = c('#ED6262','#F3ADAD','#9C2A2A')

# Create the boxplot
shift_boxplot <- ggplot(filtered_data, aes(x = condition, y = prox_to_dist_shift, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = top2_colours,labels = c("ADMAi", "SDMAi", "SDMAi+ADMAi")) +
  ylab('Proximal to Distal Usage Shift') +
  scale_x_discrete(labels=c("ADMAi" = "ADMAi", "SDMAi" = "SDMAi","SDMAi_ADMAi" = "SDMAi+ADMAi")) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.8) +
  theme_classic() +
  ylim(-0.48, 0.48) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 20)
  ) +
  labs(fill = "Condition")

# Save the plot to a file
ggsave("proximal_to_distal_usage_shift_boxplot.jpg", plot = shift_boxplot, width = 10, height = 8, dpi = 300)
