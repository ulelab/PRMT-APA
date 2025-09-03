# Load required libraries
library(tidyverse)   # For data manipulation and visualization (includes dplyr, ggplot2)
library(data.table)  # For efficient data manipulation
library(RColorBrewer)# For color palettes in plots
library(ggh4x)       # Extensions for ggplot2
library(stringr)     # For string manipulation

# ----------------------------
# Set Working Directories and File Paths
# ----------------------------

# Update these paths to match your data locations
polyA_csv_dir <- "/path/to/your/CSV_files"
anno_bed_dir <- "/path/to/your/annotated_BED_files"
UTR_info_path <- "/path/to/your/UTR_regions.bed"
splice_site_info_path <- "/path/to/your/splice_regions.bed"
output_dir <- "/path/to/your/output_directory"

# For demonstration purposes, replace the above paths with your actual paths
# For example:
polyA_csv_dir <- "/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/CSV_files"
anno_bed_dir <- "/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/TUTR_vs_ALE"
UTR_info_path <- "/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/no_dedup/two_step_test/TUTR_vs_ALE/UTR_regions_gencode_v45.bed"
splice_site_info_path <- "/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/no_dedup/two_step_test/TUTR_vs_ALE/splice_regions_gencode_v45.bed"
output_dir <- "/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/TUTR_vs_ALE"

# Ensure the output directory exists
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the experimental conditions
conditions <- c("siCFIM25_DMAi", "siCFIM25_DMSO")  # Update this vector as per your conditions

# ----------------------------
# Read and Prepare PolyA Site Data
# ----------------------------

# Set working directory to where your PolyA site CSV files are located
setwd(polyA_csv_dir)

# List all CSV files containing 'all' in the current working directory
polyA_files <- list.files(pattern = "all.*\\.csv$", full.names = TRUE)

# Read each file into a list of data frames
polyA_data_list <- lapply(polyA_files, read.csv)

# Tidy up 'ensembl_ID' column by removing version numbers or any suffix after the dot
polyA_data_list <- lapply(polyA_data_list, function(df){
  df$ensembl_ID <- sub("\\..*$", "", df$ensembl_ID)
  return(df)
})

# Calculate mean usage of each PolyA site across control and treatment conditions
polyA_data_list <- lapply(polyA_data_list, function(df){
  df$mean_usage <- (df$pA_usage_control + df$pA_usage_treatment) / 2 
  return(df)
})

# Function to process PolyA site data and select top 2 sites per gene
process_gene_groups <- function(df) {
  # Add a flag for significant PolyA sites based on adjusted p-value and change in usage
  df <- df %>%
    mutate(is_significant = ifelse(!is.na(twostep_feature_padj) & 
                                     twostep_feature_padj < 0.05 & 
                                     abs(change_in_usage) > 0.1, TRUE, FALSE))
  
  # Count the number of significant PolyA sites per gene
  df <- df %>%
    group_by(gene_name) %>%
    mutate(sig_PAS_count = sum(is_significant, na.rm = TRUE)) %>%
    ungroup()
  
  # For genes with at least two significant PAS
  df_sig <- df %>%
    dplyr::filter(sig_PAS_count >= 2 & is_significant) %>%  # Keep only significant PAS
    group_by(gene_name) %>%
    arrange(desc(abs(change_in_usage))) %>%  # Order by absolute change in usage
    slice_head(n = 2) %>%  # Select top 2 PAS per gene
    ungroup()
  
  # For genes with fewer than two significant PAS
  df_non_sig <- df %>%
    dplyr::filter(sig_PAS_count < 2) %>%
    group_by(gene_name) %>%
    arrange(desc(mean_usage)) %>%  # Order by mean usage
    slice_head(n = 2) %>%  # Select top 2 PAS per gene
    ungroup()
  
  # Combine the two data frames
  df_result <- bind_rows(df_sig, df_non_sig) %>% dplyr::select(-is_significant)
  
  return(df_result)
}

# Apply the function to each data frame in the list
polyA_top2_list <- lapply(polyA_data_list, process_gene_groups)

# ----------------------------
# Read and Prepare Annotated BED Files
# ----------------------------

# Set working directory to where your annotated BED files are located
setwd(anno_bed_dir)

# Get a list of all BED files in the directory
bed_files <- list.files(pattern = '\\.bed$', full.names = TRUE)

# Read each BED file into a list of data frames
anno_data_list <- lapply(bed_files, function(file) {
  read.delim(file, header = FALSE, stringsAsFactors = FALSE,
             col.names = c("chr", "start", "end", "col4", "gene_name", "strand",
                           "gtf_chr", "database", "region", "gtf_start",
                           "gtf_end", "col12", "gtf_strand", "col14", "gtf_annotation"))
})

# Separate 'gtf_annotation' column to extract 'ensembl_ID' and 'extra_annotation'
anno_data_list <- lapply(anno_data_list, function(df) {
  df %>%
    separate(gtf_annotation, into = c("ensembl_ID", "extra_annotation"), sep = ";", extra = "merge", fill = "right") %>%
    mutate(
      ensembl_ID = gsub("^gene_id ", "", ensembl_ID),  # Remove 'gene_id ' prefix
      ensembl_ID = sub("\\..*", "", ensembl_ID)        # Remove version numbers
    )
})

# Ensure that 'polyA_top2_list' and 'anno_data_list' are in the same order
# Join the processed PolyA data with the annotated sites
anno_data_list <- Map(function(df1, df2) {
  dplyr::left_join(df1, df2, by = c('chr', 'start', 'end', 'strand', 'gene_name', 'ensembl_ID')) 
}, polyA_top2_list, anno_data_list)

# Further separate 'extra_annotation' to extract 'transcript_annotation' and 'further_annotation'
anno_data_list <- lapply(anno_data_list, function(df) {
  df %>%
    mutate(
      transcript_annotation = if_else(
        grepl("transcript_id", extra_annotation),
        sub(";.*", "", extra_annotation),
        extra_annotation
      ),
      further_annotation = if_else(
        grepl("transcript_id", extra_annotation),
        sub(".*?;", "", extra_annotation),
        NA_character_
      )
    ) %>%
    dplyr::select(-extra_annotation)
})

# If UTR annotation is present in each 'transcript_annotation' group, only keep rows with 'UTR' in 'region'
anno_data_UTR_filtered <- lapply(anno_data_list, function(df) {
  df %>%
    group_by(chr, start, end, transcript_annotation) %>%
    dplyr::filter(!any(grepl("UTR", region)) | grepl("UTR", region)) %>%
    ungroup()
})

# Remove 'gene' annotations if there are other annotations for that PolyA site
anno_data_UTR_filtered <- lapply(anno_data_UTR_filtered, function(df) {
  df %>%
    group_by(chr, start, end) %>%
    mutate(has_extra_anno = any(region != 'gene')) %>%
    dplyr::filter(!(region == 'gene' & has_extra_anno)) %>%
    dplyr::select(-has_extra_anno) %>%
    ungroup()
})

# Remove 'transcript' annotations if there are other annotations for that PolyA site
anno_data_UTR_filtered <- lapply(anno_data_UTR_filtered, function(df) {
  df %>%
    group_by(chr, start, end) %>%
    mutate(has_extra_anno = any(region != 'transcript')) %>%
    dplyr::filter(!(region == 'transcript' & has_extra_anno)) %>%
    dplyr::select(-has_extra_anno) %>%
    ungroup()
})

# Identify genes that are associated with more than one PolyA site
filtered_genes_list <- lapply(anno_data_UTR_filtered, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    summarize(
      distinct_site_count = n_distinct(paste(chr, start, end, sep = "_"))
    ) %>%
    dplyr::filter(distinct_site_count > 1) %>%
    dplyr::select(ensembl_ID) %>%
    ungroup()
})

# Keep only the genes that are associated with multiple PolyA sites
anno_data_single_site_filtered <- Map(function(df1, df2) {
  semi_join(df1, df2, by = "ensembl_ID")
}, anno_data_UTR_filtered, filtered_genes_list)

# Find sites that were completely lost due to the filtering and restore them
lost_sites_list <- Map(function(df1, df2) {
  anti_join(df1, df2, by = c("chr", "start", "end"))
}, anno_data_UTR_filtered, anno_data_single_site_filtered)

# Combine the filtered data frames with the lost sites
anno_data_UTR_filtered <- Map(rbind, anno_data_single_site_filtered, lost_sites_list)

# Read in UTR location file
UTR_info <- read.delim(UTR_info_path, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("UTR_chr", "UTR_start", "UTR_end", "UTR_strand", "ensembl_ID", "UTR_class", "transcript_annotation"))

# Keep only 3' UTR information
UTR3_info <- UTR_info %>% 
  dplyr::filter(UTR_class == '3UTR')

# Prepare UTR3_info for merging
UTR3_info$UTR_chr <- paste0("chr", UTR3_info$UTR_chr)
UTR3_info$ensembl_ID <- sub("\\..*$", "", UTR3_info$ensembl_ID)
UTR3_info$transcript_annotation <- paste(" transcript_id ", UTR3_info$transcript_annotation, sep = "")

# Merge annotated sites with UTR start and end information
anno_data_with_UTR <- lapply(anno_data_UTR_filtered, function(df) {
  dplyr::left_join(df, UTR3_info, by = c("ensembl_ID", "transcript_annotation"))
})

# Split data frames into positive and negative strands
anno_data_pos <- lapply(anno_data_with_UTR, function(df) {
  df %>% dplyr::filter(strand == '+')
})
anno_data_neg <- lapply(anno_data_with_UTR, function(df) {
  df %>% dplyr::filter(strand == '-')
})

# For positive strand, keep annotations where transcript end is near PolyA site
anno_data_pos_filtered <- lapply(anno_data_pos, function(df) {
  df %>%
    mutate(
      keep_row = ifelse(is.na(UTR_end),
                        gtf_end > end - 100 & gtf_end < end + 30,
                        UTR_end > end - 100 & UTR_end < end + 30)
    ) %>%
    dplyr::filter(keep_row) %>%
    dplyr::select(-keep_row)
})

# For negative strand
anno_data_neg_filtered <- lapply(anno_data_neg, function(df) {
  df %>%
    mutate(
      keep_row = ifelse(is.na(UTR_end),
                        gtf_start < start + 100 & gtf_start > start - 30,
                        UTR_start < start + 100 & UTR_start > start - 30)
    ) %>%
    dplyr::filter(keep_row) %>%
    dplyr::select(-keep_row)
})

# Function to find lost sites
find_lost_sites <- function(original_df, filtered_df) {
  anti_join(original_df, filtered_df, by = c("chr", "start", "end"))
}

# Find lost sites for positive and negative strands
lost_sites_pos <- Map(find_lost_sites, anno_data_pos, anno_data_pos_filtered)
lost_sites_neg <- Map(find_lost_sites, anno_data_neg, anno_data_neg_filtered)

# For lost sites on positive strand, keep the transcript with end closest to the PolyA site
lost_sites_pos_filtered <- lapply(lost_sites_pos, function(df) {
  df %>%
    mutate(further_annotation = coalesce(further_annotation, 'None')) %>%
    group_by(chr, start, end, ensembl_ID) %>%
    dplyr::filter(
      if (any(str_detect(further_annotation, 'transcript_support_level 1'))) {
        str_detect(further_annotation, 'transcript_support_level 1') | str_detect(further_annotation, 'transcript_support_level 2')
      } else if (any(str_detect(further_annotation, 'transcript_support_level 2'))) {
        str_detect(further_annotation, 'transcript_support_level 2')
      } else if (any(str_detect(further_annotation, 'transcript_support_level 3'))) {
        str_detect(further_annotation, 'transcript_support_level 3')
      } else {
        TRUE
      }
    ) %>%
    ungroup() %>%
    mutate(
      dist_to_end = ifelse(is.na(UTR_end), abs(gtf_end - end), abs(UTR_end - end))
    ) %>%
    group_by(chr, start, end, ensembl_ID) %>%
    dplyr::filter(dist_to_end == min(dist_to_end)) %>%
    dplyr::select(-dist_to_end) %>%
    ungroup()
})

# For lost sites on negative strand
lost_sites_neg_filtered <- lapply(lost_sites_neg, function(df) {
  df %>%
    mutate(further_annotation = coalesce(further_annotation, 'None')) %>%
    group_by(chr, start, end, ensembl_ID) %>%
    dplyr::filter(
      if (any(str_detect(further_annotation, 'transcript_support_level 1'))) {
        str_detect(further_annotation, 'transcript_support_level 1') | str_detect(further_annotation, 'transcript_support_level 2')
      } else if (any(str_detect(further_annotation, 'transcript_support_level 2'))) {
        str_detect(further_annotation, 'transcript_support_level 2')
      } else if (any(str_detect(further_annotation, 'transcript_support_level 3'))) {
        str_detect(further_annotation, 'transcript_support_level 3')
      } else {
        TRUE
      }
    ) %>%
    ungroup() %>%
    mutate(
      dist_to_end = ifelse(is.na(UTR_start), abs(gtf_start - start), abs(UTR_start - start))
    ) %>%
    group_by(chr, start, end, ensembl_ID) %>%
    dplyr::filter(dist_to_end == min(dist_to_end)) %>%
    dplyr::select(-dist_to_end) %>%
    ungroup()
})

# Combine the filtered data frames with the lost sites
anno_data_pos_combined <- Map(rbind, anno_data_pos_filtered, lost_sites_pos_filtered)
anno_data_neg_combined <- Map(rbind, anno_data_neg_filtered, lost_sites_neg_filtered)

# Combine positive and negative strands
anno_data_combined <- Map(rbind, anno_data_pos_combined, anno_data_neg_combined)

# Filter coordinates by transcript support level (TSL) again
anno_data_TSL_filtered <- lapply(anno_data_combined, function(df) {
  df %>%
    mutate(further_annotation = coalesce(further_annotation, 'None')) %>%
    group_by(chr, start, end, ensembl_ID) %>%
    dplyr::filter(
      if (any(str_detect(further_annotation, 'transcript_support_level 1'))) {
        str_detect(further_annotation, 'transcript_support_level 1') | str_detect(further_annotation, 'transcript_support_level 2')
      } else if (any(str_detect(further_annotation, 'transcript_support_level 2'))) {
        str_detect(further_annotation, 'transcript_support_level 2')
      } else if (any(str_detect(further_annotation, 'transcript_support_level 3'))) {
        str_detect(further_annotation, 'transcript_support_level 3')
      } else {
        TRUE
      }
    ) %>%
    ungroup()
})

# Read in splice site location file
splice_site_info <- read.delim(splice_site_info_path, header = FALSE, stringsAsFactors = FALSE,
                               col.names = c("splice_chr", "splice_start", "splice_end", "coordinates", "splice_score", "splice_strand", "splice_site_type", "transcript_annotation", "ensembl_ID"))

# Keep only final splice site for each transcript
final_splice_site_info <- splice_site_info %>%
  group_by(transcript_annotation) %>%
  dplyr::filter(ifelse(splice_strand == '+', splice_start == max(splice_start), splice_end == min(splice_end))) %>%
  ungroup()

# Prepare splice site info for merging
final_splice_site_info$splice_chr <- paste0("chr", final_splice_site_info$splice_chr)
final_splice_site_info$ensembl_ID <- sub("\\..*$", "", final_splice_site_info$ensembl_ID)
final_splice_site_info$transcript_annotation <- paste(" transcript_id ", final_splice_site_info$transcript_annotation, sep = "")

# Merge annotated sites with splice site information
sites_with_splice_info_list <- lapply(anno_data_TSL_filtered, function(df) {
  dplyr::left_join(df, final_splice_site_info, by = c("ensembl_ID", "transcript_annotation"))
})

# Change NA 'splice_start' values to 'None'
sites_with_splice_info_list <- lapply(sites_with_splice_info_list, function(df) {
  df$splice_start <- ifelse(is.na(df$splice_start), 'None', df$splice_start)
  return(df)
})

# Add 'UTR_type' column based on criteria
sites_with_UTR_type_list <- lapply(sites_with_splice_info_list, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    mutate(
      UTR_type = case_when(
        n_distinct(transcript_annotation) == 1 ~ 'TUTR',
        all(is.na(splice_start)) ~ 'TUTR',
        TRUE ~ 'UNKNOWN'
      )
    ) %>%
    ungroup()
})

# Separate 'TUTR' sites
TUTR_sites_list <- lapply(sites_with_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'TUTR')
})

# Identify 'single_sig_pA' genes from the 'TUTR' sites
TUTR_sites_list <- lapply(TUTR_sites_list, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    mutate(
      UTR_type = case_when(
        n_distinct(start) == 1 ~ 'single_sig_pA',
        TRUE ~ 'TUTR'
      )
    ) %>%
    ungroup()
})

# Process 'UNKNOWN' sites
unknown_sites_list <- lapply(sites_with_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'UNKNOWN')
})

# If there is only one terminal splice site per gene, classify as 'TUTR'
sites_with_UTR_type_list <- lapply(unknown_sites_list, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    mutate(
      UTR_type = case_when(
        n_distinct(splice_start) == 1 ~ 'TUTR',
        TRUE ~ 'UNKNOWN'
      )
    ) %>%
    ungroup()
})

# Add newly identified 'TUTR' sites
more_TUTR_sites_list <- lapply(sites_with_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'TUTR')
})

# Check for 'single_pA' again
more_TUTR_sites_list <- lapply(more_TUTR_sites_list, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    mutate(
      UTR_type = case_when(
        n_distinct(start) == 1 ~ 'single_pA',
        TRUE ~ 'TUTR'
      )
    ) %>%
    ungroup()
})

# Combine 'TUTR' sites
TUTR_sites_list <- Map(rbind, TUTR_sites_list, more_TUTR_sites_list)

# Identify 'iAPA' sites within 'TUTR' sites
TUTR_sites_list <- lapply(TUTR_sites_list, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    mutate(
      UTR_type = case_when(
        strand == "+" & splice_start >= min(start, na.rm = TRUE) & splice_start <= max(start, na.rm = TRUE) ~ "iAPA",
        strand == "-" & splice_end >= min(end, na.rm = TRUE) & splice_end <= max(end, na.rm = TRUE) ~ "iAPA",
        TRUE ~ UTR_type
      )
    ) %>%
    ungroup()
})

# Separate 'iAPA' sites
iAPA_sites_list <- lapply(TUTR_sites_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'iAPA')
})

# Update 'TUTR' sites
TUTR_sites_list <- lapply(TUTR_sites_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'TUTR')
})

# Process remaining 'UNKNOWN' sites
remaining_unknown_sites_list <- lapply(sites_with_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'UNKNOWN')
})

# Classify as 'ALE' or 'MIXED' based on criteria
sites_with_ALE_UTR_type_list <- lapply(remaining_unknown_sites_list, function(df) {
  df %>%
    group_by(ensembl_ID, splice_start) %>%
    mutate(num_unique_starts = n_distinct(start)) %>%
    ungroup() %>%
    group_by(ensembl_ID) %>%
    mutate(
      has_multiple_starts = any(num_unique_starts > 1),
      UTR_type = ifelse(has_multiple_starts, "MIXED", "ALE")
    ) %>%
    dplyr::select(-has_multiple_starts, -num_unique_starts) %>%
    ungroup()
})

# Separate 'ALE' sites
ALE_sites_list <- lapply(sites_with_ALE_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'ALE')
})

# Separate 'MIXED' sites
MIXED_sites_list <- lapply(sites_with_ALE_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'MIXED')
})

# Combine all sites into one list
all_sites_with_UTR_type <- Map(rbind, TUTR_sites_list, sites_with_ALE_UTR_type_list)
all_sites_with_UTR_type <- Map(rbind, all_sites_with_UTR_type, iAPA_sites_list)

# Keep relevant columns and remove duplicate rows
all_sites_with_UTR_type <- lapply(all_sites_with_UTR_type, function(df) {
  df <- df[, c(1:18, 42), drop = FALSE] %>% distinct()
  return(df)
})

# ----------------------------
# Filter Significant and Control Genes
# ----------------------------

# Filter to isolate significant genes (genes with at least two significant PAS)
sig_sites_with_UTR_type <- lapply(all_sites_with_UTR_type, function(df) {
  df %>% dplyr::filter(sig_PAS_count >= 2)
})

# Filter to isolate control genes (genes with stable PAS usage)
ctrl_sites_with_UTR_type <- lapply(all_sites_with_UTR_type, function(df) {
  df %>%
    dplyr::filter(mean_usage > 0.075 & abs(change_in_usage) < 0.075) %>%
    group_by(gene_name) %>%
    dplyr::filter(n() >= 2) %>%
    ungroup()
})

# ----------------------------
# Counting Genes per UTR Type and Plotting
# ----------------------------

# Function to count genes for each 'UTR_type'
count_genes_by_UTR_type <- function(df_list, target_UTR_type) {
  lapply(df_list, function(df) {
    df %>%
      group_by(gene_name) %>%
      dplyr::filter(any(UTR_type == target_UTR_type)) %>%
      ungroup() %>%
      summarise(num_genes = n_distinct(gene_name))
  })
}

# Function to safely extract 'num_genes' from each dataframe
extract_num_genes <- function(df) {
  if (nrow(df) == 0) {
    return(0)
  } else {
    return(df$num_genes)
  }
}

# Count number of genes belonging to each 'UTR_type' for significant genes
sig_TUTR_genes_count <- count_genes_by_UTR_type(sig_sites_with_UTR_type, "TUTR")
sig_ALE_genes_count <- count_genes_by_UTR_type(sig_sites_with_UTR_type, "ALE")
sig_MIXED_genes_count <- count_genes_by_UTR_type(sig_sites_with_UTR_type, "MIXED")
sig_iAPA_genes_count <- count_genes_by_UTR_type(sig_sites_with_UTR_type, "iAPA")

# Create a data frame with counts for each 'UTR_type' and each condition
UTR_gene_tally_list <- lapply(1:length(sig_TUTR_genes_count), function(i) {
  data.frame(
    UTR_type = c('TUTR', 'ALE', 'MIXED', 'iAPA'),
    total = c(
      extract_num_genes(sig_TUTR_genes_count[[i]]),
      extract_num_genes(sig_ALE_genes_count[[i]]),
      extract_num_genes(sig_MIXED_genes_count[[i]]),
      extract_num_genes(sig_iAPA_genes_count[[i]])
    ),
    condition = conditions[i]
  )
})

# Combine all tallies into one data frame for plotting
combined_tallies <- do.call(rbind, UTR_gene_tally_list)

# Convert 'UTR_type' to factor with specified levels for consistent plotting
combined_tallies$UTR_type <- factor(combined_tallies$UTR_type, levels = c('TUTR', 'ALE', 'MIXED', 'iAPA'))

# Define custom labels for the strip titles
strip_labels <- c(
  "siCFIM25_DMAi"      = "siCFIM25 + DMAi",
  "siCFIM25_DMSO"      = "siCFIM25 + DMSO"
  # Add more mappings as needed
)

# Plot the bar chart for significant genes
sig_UTR_class_barchart <- ggplot(combined_tallies, aes(x = UTR_type, y = total, fill = UTR_type)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.7) +
  facet_wrap(~condition, scales = "free_x", labeller = labeller(condition = strip_labels)) +
  labs(
    x = 'UTR Type',
    y = 'Number of Genes'
  ) +
  theme_classic() +
  scale_fill_brewer(palette = 'Set2') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgrey", color = "black", size = 1)
  )

# Save the bar chart plot to a file
ggsave(filename = file.path(output_dir, "sig_UTR_class_barchart.jpg"), 
       plot = sig_UTR_class_barchart, width = 10, height = 8, dpi = 300)

# ----------------------------
# Repeat Counting and Plotting for Control Genes
# ----------------------------

# Count number of genes belonging to each 'UTR_type' for control genes
ctrl_TUTR_genes_count <- count_genes_by_UTR_type(ctrl_sites_with_UTR_type, "TUTR")
ctrl_ALE_genes_count <- count_genes_by_UTR_type(ctrl_sites_with_UTR_type, "ALE")
ctrl_MIXED_genes_count <- count_genes_by_UTR_type(ctrl_sites_with_UTR_type, "MIXED")
ctrl_iAPA_genes_count <- count_genes_by_UTR_type(ctrl_sites_with_UTR_type, "iAPA")

# Create a data frame with counts for each 'UTR_type' and each condition
UTR_gene_tally_ctrl_list <- lapply(1:length(ctrl_TUTR_genes_count), function(i) {
  data.frame(
    UTR_type = c('TUTR', 'ALE', 'MIXED', 'iAPA'),
    total = c(
      extract_num_genes(ctrl_TUTR_genes_count[[i]]),
      extract_num_genes(ctrl_ALE_genes_count[[i]]),
      extract_num_genes(ctrl_MIXED_genes_count[[i]]),
      extract_num_genes(ctrl_iAPA_genes_count[[i]])
    ),
    condition = conditions[i]
  )
})

# Combine all tallies into one data frame for plotting
combined_tallies_ctrl <- do.call(rbind, UTR_gene_tally_ctrl_list)

# Convert 'UTR_type' to factor with specified levels for consistent plotting
combined_tallies_ctrl$UTR_type <- factor(combined_tallies_ctrl$UTR_type, levels = c('TUTR', 'ALE', 'MIXED', 'iAPA'))

# Plot the bar chart for control genes
ctrl_UTR_class_barchart <- ggplot(combined_tallies_ctrl, aes(x = UTR_type, y = total, fill = UTR_type)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.7) +
  facet_wrap(~condition, scales = "free_x", labeller = labeller(condition = strip_labels)) +
  labs(
    x = 'UTR Type',
    y = 'Number of Genes',
  ) +
  theme_classic() +
  scale_fill_brewer(palette = 'Set2') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgrey", color = "black", size = 1)
  )

# Save the bar chart plot to a file
ggsave(filename = file.path(output_dir, "ctrl_UTR_class_barchart.jpg"), 
       plot = ctrl_UTR_class_barchart, width = 10, height = 8, dpi = 300)

# ----------------------------
# Jitter Plot for Proximal to Distal Usage Shift
# ----------------------------

# Define functions to identify positive and negative strands
identify_positive_strands <- function(df) {
  df %>% dplyr::filter(strand == "+")
}

identify_negative_strands <- function(df) {
  df %>% dplyr::filter(strand == "-")
}

# Define functions to identify distal and proximal PolyA sites
identify_distal_pAs <- function(pos_df, neg_df) {
  pos_distal <- pos_df %>% group_by(gene_name) %>% dplyr::slice_max(order_by = start, n = 1, with_ties = FALSE)
  neg_distal <- neg_df %>% group_by(gene_name) %>% dplyr::slice_min(order_by = start, n = 1, with_ties = FALSE)
  bind_rows(pos_distal, neg_distal)
}

identify_proximal_pAs <- function(pos_df, neg_df) {
  pos_proximal <- pos_df %>% group_by(gene_name) %>% dplyr::slice_min(order_by = start, n = 1, with_ties = FALSE)
  neg_proximal <- neg_df %>% group_by(gene_name) %>% dplyr::slice_max(order_by = start, n = 1, with_ties = FALSE)
  bind_rows(pos_proximal, neg_proximal)
}

# **Step 1:** Identify positive and negative strand transcripts for significant genes
sig_pos_sites_with_UTR_type_anno_list <- lapply(sig_sites_with_UTR_type, identify_positive_strands)
sig_neg_sites_with_UTR_type_anno_list <- lapply(sig_sites_with_UTR_type, identify_negative_strands)
ctrl_pos_sites_with_UTR_type_anno_list <- lapply(ctrl_sites_with_UTR_type, identify_positive_strands)
ctrl_neg_sites_with_UTR_type_anno_list <- lapply(ctrl_sites_with_UTR_type, identify_negative_strands)

# **Step 2:** Identify distal and proximal PolyA sites
sig_anno_distal_pAs_list <- Map(identify_distal_pAs, sig_pos_sites_with_UTR_type_anno_list, sig_neg_sites_with_UTR_type_anno_list)
sig_anno_proximal_pAs_list <- Map(identify_proximal_pAs, sig_pos_sites_with_UTR_type_anno_list, sig_neg_sites_with_UTR_type_anno_list)
ctrl_anno_distal_pAs_list <- Map(identify_distal_pAs, ctrl_pos_sites_with_UTR_type_anno_list, ctrl_neg_sites_with_UTR_type_anno_list)
ctrl_anno_proximal_pAs_list <- Map(identify_proximal_pAs, ctrl_pos_sites_with_UTR_type_anno_list, ctrl_neg_sites_with_UTR_type_anno_list)

# **Step 3:** Combine distal and proximal PolyA sites with 'pA_type' column
sig_sites_with_UTR_type_anno_with_position_list <- Map(function(dist_df, prox_df) {
  dist_df$pA_type <- "distal pA"
  prox_df$pA_type <- "proximal pA"
  bind_rows(dist_df, prox_df)
}, sig_anno_distal_pAs_list, sig_anno_proximal_pAs_list)

ctrl_sites_with_UTR_type_anno_with_position_list <- Map(function(dist_df, prox_df) {
  dist_df$pA_type <- "distal pA"
  prox_df$pA_type <- "proximal pA"
  bind_rows(dist_df, prox_df)
}, ctrl_anno_distal_pAs_list, ctrl_anno_proximal_pAs_list)

# **Step 4:** Extract specific UTR types for significant genes 
sig_TUTR_genes_list <- lapply(sig_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "TUTR")
  df_filtered$UTR_type <- "TUTR"  # Ensure consistency
  return(df_filtered)
})

sig_ALE_genes_list <- lapply(sig_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "ALE")
  df_filtered$UTR_type <- "ALE"
  return(df_filtered)
})

sig_MIXED_genes_list <- lapply(sig_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "MIXED")
  df_filtered$UTR_type <- "MIXED"
  return(df_filtered)
})

sig_iAPA_genes_list <- lapply(sig_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "iAPA")
  df_filtered$UTR_type <- "iAPA"
  return(df_filtered)
})

# **Step 5:** Add 'condition' information to each dataframe
for (i in seq_along(sig_sites_with_UTR_type_anno_with_position_list)) {
  sig_sites_with_UTR_type_anno_with_position_list[[i]]$condition <- conditions[i]
}

# **Step 6:** Combine all UTR types into one data frame per condition
sig_prox_dist_shift_list <- list()

for (i in seq_along(sig_sites_with_UTR_type_anno_with_position_list)) {
  # Extract data frames
  all_df <- sig_sites_with_UTR_type_anno_with_position_list[[i]]  # all UTR types
  TUTR_df <- sig_TUTR_genes_list[[i]]                             # TUTR
  ALE_df <- sig_ALE_genes_list[[i]]                               # ALE
  MIXED_df <- sig_MIXED_genes_list[[i]]                           # MIXED
  iAPA_df <- sig_iAPA_genes_list[[i]]                             # iAPA
  
  # Combine all UTR types into one data frame
  combined_df <- bind_rows(all_df, TUTR_df, ALE_df, MIXED_df, iAPA_df)
  
  # Identify distal and proximal PolyA sites
  pos_df <- combined_df %>% dplyr::filter(strand == "+")
  neg_df <- combined_df %>% dplyr::filter(strand == "-")
  
  distal_pAs <- identify_distal_pAs(pos_df, neg_df)
  proximal_pAs <- identify_proximal_pAs(pos_df, neg_df)
  
  # Add 'pA_type' and combine
  distal_pAs$pA_type <- "distal pA"
  proximal_pAs$pA_type <- "proximal pA"
  prox_dist_combined <- bind_rows(distal_pAs, proximal_pAs)
  
  # Calculate 'prox_to_dist_shift' per gene, UTR_type, and condition
  prox_dist_shift <- prox_dist_combined %>%
    group_by(gene_name, UTR_type, condition) %>%
    mutate(prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] - change_in_usage[pA_type == "proximal pA"]) %>%
    ungroup()
  
  # Remove outlier genes based on 'prox_to_dist_shift'
  prox_dist_shift_filtered <- prox_dist_shift %>% 
    dplyr::filter(abs(prox_to_dist_shift) > 0.2)
  
  # Append to the list
  sig_prox_dist_shift_list[[i]] <- prox_dist_shift_filtered
}

# **Step 7:** Combine all processed data into one data frame
combined_prox_dist_shift <- bind_rows(sig_prox_dist_shift_list)

# Duplicate the data frame and change UTR_type to 'all' in the duplicate
all_utr_df <- combined_prox_dist_shift %>%
  mutate(UTR_type = 'all')

# Combine the original data frame with the modified duplicate
combined_prox_dist_shift <- bind_rows(combined_prox_dist_shift, all_utr_df)

# **Step 8:** Ensure 'UTR_type' and 'pA_type' are factors with specified levels
combined_prox_dist_shift$UTR_type <- factor(combined_prox_dist_shift$UTR_type, levels = c('all','TUTR', 'ALE', 'MIXED', 'iAPA'))

# Get the first four colors from the Set2 palette
set2_colors <- brewer.pal(8, "Set2")[1:4]

# Construct the color vector
color_vector <- c(
  '#ED6262',          # First color
  set2_colors,        # First four elements from Set2
  '#F3ADAD',          # Middle color
  set2_colors,        # Repeat first four elements from Set2
  '#9C2A2A',          # Last color
  set2_colors         # Repeat first four elements from Set2 again
)

# Create a unique identifier for each combination of UTR_type and condition
combined_prox_dist_shift$group_id <- with(combined_prox_dist_shift, interaction(UTR_type, condition))
# Check how many unique combinations there are
unique_combinations <- unique(combined_prox_dist_shift$group_id)
# Create a mapping from combinations to colors
color_mapping <- setNames(color_vector, unique_combinations)
# Map the colors to your dataframe
combined_prox_dist_shift$colors <- color_mapping[combined_prox_dist_shift$group_id]

#jitter plot
jitter_plot = ggplot(combined_prox_dist_shift, aes(x = UTR_type, y = prox_to_dist_shift, fill = colors)) +
  geom_violin(trim = FALSE, alpha = 0.4) + # Add violin plot to show density
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.3), 
              size = 1.5, aes(color = colors), alpha = 0.6) + # Add jitter points for visibility
  stat_summary(fun = median, geom = "point", 
               position = position_dodge(width = 0.3), 
               size = 6, shape = 22, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.8) +
  facet_wrap(~condition, scales = "free_x", labeller = labeller(condition = strip_labels)) +
  scale_fill_identity() + # Use the actual fill colors provided in the 'colors' attribute
  scale_color_identity() + # Use the actual colors provided in the 'colors' attribute
  labs(
    x = 'UTR Type',
    y = "Proximal to Distal Usage Shift",
    color = 'UTR Type',
    fill = 'UTR Type'
  ) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 28),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgrey", color = "black", size = 1)
  ) +
  ylim(-1.5, 1.5)

# Save the jitter plot to a file
ggsave(filename = file.path(output_dir, "proximal_to_distal_usage_shift_TUTR_vs_ALE_jitter_plot.jpg"), 
       plot = jitter_plot, width = 10, height = 8, dpi = 300)

# ----------------------------
# Write Significant and Control Data Frames to CSV Files
# ----------------------------

# Prepare significant data frames for writing
for (i in seq_along(sig_sites_with_UTR_type)) {
  condition <- conditions[i]
  
  # Create lists of data frames for each condition
  sig_list <- list(
    sig_sites_with_UTR_type[[i]],  # all
    sig_TUTR_genes_list[[i]],      # TUTR
    sig_ALE_genes_list[[i]],       # ALE
    sig_MIXED_genes_list[[i]],     # MIXED
    sig_iAPA_genes_list[[i]]       # iAPA
  )
  
  # Define the directory name
  output_dir <- "sig_CSV/"
  
  # Create the directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Define file names as in the previous function
  file_names <- c(
    paste0("sig_", condition, "_all.csv"),
    paste0("sig_", condition, "_TUTR.csv"),
    paste0("sig_", condition, "_ALE.csv"),
    paste0("sig_", condition, "_MIXED.csv"),
    paste0("sig_", condition, "_iAPA.csv")
  )
  
  # Write each data frame to a CSV file in the new directory
  for (j in seq_along(sig_list)) {
    file_path <- file.path(output_dir, file_names[j])  # Combine directory and file name
    write.csv(sig_list[[j]], file = file_path, row.names = FALSE)
  }
}

#same for controls
ctrl_TUTR_genes_list <- lapply(ctrl_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "TUTR")
  df_filtered$UTR_type <- "TUTR"  # Ensure consistency
  return(df_filtered)
})

ctrl_ALE_genes_list <- lapply(ctrl_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "ALE")
  df_filtered$UTR_type <- "ALE"
  return(df_filtered)
})

ctrl_MIXED_genes_list <- lapply(ctrl_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "MIXED")
  df_filtered$UTR_type <- "MIXED"
  return(df_filtered)
})

ctrl_iAPA_genes_list <- lapply(ctrl_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "iAPA")
  df_filtered$UTR_type <- "iAPA"
  return(df_filtered)
})

# Prepare control data frames for writing
for (i in seq_along(ctrl_sites_with_UTR_type)) {
  condition <- conditions[i]
  
  # Create lists of data frames for each condition
  ctrl_list <- list(
    ctrl_sites_with_UTR_type[[i]],  # all
    ctrl_TUTR_genes_list[[i]],      # TUTR
    ctrl_ALE_genes_list[[i]],       # ALE
    ctrl_MIXED_genes_list[[i]],     # MIXED
    ctrl_iAPA_genes_list[[i]]       # iAPA
  )
  
  # Define the directory name
  output_dir <- "ctrl_CSV/"
  
  # Create the directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Define file names as in the previous function
  file_names <- c(
    paste0("ctrl_", condition, "_all.csv"),
    paste0("ctrl_", condition, "_TUTR.csv"),
    paste0("ctrl_", condition, "_ALE.csv"),
    paste0("ctrl_", condition, "_MIXED.csv"),
    paste0("ctrl_", condition, "_iAPA.csv")
  )
  
  # Write each data frame to a CSV file in the new directory
  for (j in seq_along(ctrl_list)) {
    file_path <- file.path(output_dir, file_names[j])  # Combine directory and file name
    write.csv(ctrl_list[[j]], file = file_path, row.names = FALSE)
  }
}
