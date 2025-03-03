library(tidyverse)
library(data.table)
library(dunn.test)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(ggpattern)
library(GenomicFeatures)


##read in DMAi pA site CSVs and merge with DESEQ2 tables to obtain transcript expression data (BaseMean)
#setwd 
setwd("/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/panel_of_cancer_lines/positional_info")

#read in mulit_sig_pAs file
csv_files = list.files(pattern = "^top2")

# List to store data frames
DMAi_data_frames = list()

# Loop through each CSV file and read it
for (file in csv_files) {
  df <- read_csv(file)
  DMAi_data_frames[[file]] <- df
}

#read in DE tables
setwd("/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/DESeq2_tables/DMAi_cancer_line_panel")

#read in DE file
DE_csv_files = list.files(pattern = "csv")

# List to store data frames
DE_data_frames = list()

# Loop through each CSV file and read it
for (file in DE_csv_files) {
  df <- read_csv(file)
  DE_data_frames[[file]] <- df
}

#keep relevant columns
DE_data_frames = lapply(DE_data_frames,function(df){
  df = df[,c('baseMean','gene_name')]
  return(df)
})

#add baseMean to tables
DMAi_data_frames_with_basemean <- Map(function(x, y) {
  merge(x, y, by = "gene_name", all.x = TRUE)
}, DMAi_data_frames, DE_data_frames)

##read in non-DMAi pA site CSVs and merge with DESEQ2 tables to obtain transcript expression data (BaseMean)
#setwd
setwd("/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig4/positional_info/")

#read in mulit_sig_pAs file
csv_files = list.files(pattern = "^top2")

# List to store data frames
non_DMAi_data_frames = list()

# Loop through each CSV file and read it
for (file in csv_files) {
  df <- read_csv(file)
  non_DMAi_data_frames[[file]] <- df
}

#read in DE file
setwd("/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/DESeq2_tables/non_DMAi")
non_DMAi_DE_csv_files = list.files(pattern = "csv")

# List to store data frames
DE_data_frames = list()

# Loop through each CSV file and read it
for (file in non_DMAi_DE_csv_files) {
  df <- read_csv(file)
  DE_data_frames[[file]] <- df
}

#keep relevant columns
DE_data_frames = lapply(DE_data_frames,function(df){
  df = df[,c('baseMean','gene_name')]
  return(df)
})

#add baseMean to tables
non_DMAi_data_frames_with_basemean <- Map(function(x, y) {
  merge(x, y, by = "gene_name", all.x = TRUE)
}, non_DMAi_data_frames, DE_data_frames)

##find common DMAi reg sites that are a) present in at least 2 cell lines, b)  avg effect > 10% c) significant in at least 1, d) change by > 7.5% in all 
# Define the columns of interest and thresholds
columns_of_interest <- c("chr", "start", "end", "strand")

#filter
filter_dmai_data <- function(dmai_data_frames, columns_of_interest) {
  # Combine all data frames into one, keeping 'df_index' as a character to avoid conversion issues
  combined_df <- bind_rows(dmai_data_frames, .id = "df_index")
  
  # Filter out rows with NA in change_in_usage only
  combined_df <- combined_df %>%
    filter(!is.na(change_in_usage))
  
  # First, compute necessary statistics and add a conditional check
  results_df <- combined_df %>%
    group_by(across(all_of(columns_of_interest)))  %>%
    summarize(
      total_count = n_distinct(df_index),
      avg_change_in_usage = mean(change_in_usage, na.rm = TRUE),
      all_above_threshold = case_when(
        avg_change_in_usage > 0.1  ~ all(change_in_usage > 0.075),
        avg_change_in_usage < -0.1 ~ all(change_in_usage < -0.075),
        TRUE ~ FALSE  # Only allow groups that meet the strict conditions to pass
      ),
      below_padj_threshold_count = sum(twostep_feature_padj < 0.05, na.rm = TRUE),
      max_abs_change_in_usage = max(abs(change_in_usage), na.rm = TRUE),
      max_change_in_usage = change_in_usage[which.max(abs(change_in_usage))],
      avg_abs_change_in_usage = abs(mean(change_in_usage, na.rm = TRUE)),
      avg_basemean = mean(baseMean, na.rm = TRUE),
      avg_ctrl_usage = mean(pA_usage_control, na.rm = TRUE),
      cell_lines = toString(unique(condition)),
      gene_name = dplyr::first(gene_name),
      feature_id = dplyr::first(feature_id),
      pA_type = dplyr::first(pA_type),
      .groups = 'drop'
    ) %>%
    filter(
      total_count > 1,
      all_above_threshold,
      below_padj_threshold_count / total_count > 0,
      avg_change_in_usage > 0.1 | avg_change_in_usage < -0.1
    )
  
  return(results_df)
}

# Example usage:
filtered_dmai_results = filter_dmai_data(DMAi_data_frames_with_basemean, columns_of_interest)

## only retain pA pairs - i.e genes that have two sites passing the threshold
common_DMAi_pairs_no_non_DMAi_filter = filtered_dmai_results %>%
  distinct(chr, start, end, strand, .keep_all = TRUE) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  group_by(gene_name, pA_type) %>%  # Group by both gene_name and pA_type
  slice_max(order_by = avg_abs_change_in_usage, n = 1) %>%  # Keep only the row with the highest avg_abs_change_in_usage for each UTR_type within each gene_name
  ungroup() %>%
  group_by(gene_name) %>%
  filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA')) %>%  # Ensure each gene group has at least one 'proximal pA' and one 'distal pA'
  mutate(sign_proximal = sign(max_change_in_usage[pA_type == 'proximal pA']),  # Determine the sign of avg_change_in_usage for proximal
         sign_distal = sign(max_change_in_usage[pA_type == 'distal pA'])) %>%  # Determine the sign of avg_change_in_usage for distal
  filter(sign_proximal != sign_distal) %>%  # Keep only those gene groups where the signs are different
  dplyr::select(-sign_proximal, -sign_distal) %>%  # Optionally remove the temporary sign columns
  ungroup()

#deal with different basemeans for different sites of same gene (because some sites may have been identified in different number of lines)
common_DMAi_pairs_no_non_DMAi_filter = common_DMAi_pairs_no_non_DMAi_filter %>%
  group_by(gene_name) %>%
  mutate(avg_basemean = mean(avg_basemean, na.rm = TRUE)) %>%
  ungroup()

#write to file 
write_csv(common_DMAi_pairs_no_non_DMAi_filter,"/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig4/common_DMAi_sites.csv")

## function to add additional filter - remove sites that are altered more in non-DMAi, i.e max dmai effect is at least 33% the effect of max non-DMAi
filter_with_non_dmai <- function(dmai_results, non_dmai_data_frames, columns_of_interest) {
  # Compute the maximum absolute change_in_usage for each group across all non-DMAi data frames
  max_abs_change_in_usage_non_dmai_per_group <- bind_rows(non_dmai_data_frames) %>%
    filter(!is.na(change_in_usage)) %>%  # Ensure there are no NAs in change_in_usage
    group_by(across(all_of(columns_of_interest))) %>%
    summarize(max_abs_change_in_usage_non_dmai = max(abs(change_in_usage), na.rm = TRUE),
              .groups = 'drop')  # Drop grouping after summarization
  
  # Join the non-DMAi max values with the DMAi results
  dmai_results <- dmai_results %>%
    left_join(max_abs_change_in_usage_non_dmai_per_group, by = columns_of_interest)
  
  # Filter the results where avg_abs_change_in_usage is at least 33% of max_abs_change_in_usage_non_dmai
  filtered_results <- dmai_results %>%
    filter(avg_abs_change_in_usage >= (0.33 * max_abs_change_in_usage_non_dmai)) 
  
  return(filtered_results)
}

# Assuming 'filtered_dmai_results' is already computed from the previous step
final_results <- filter_with_non_dmai(filtered_dmai_results, non_DMAi_data_frames, columns_of_interest)

common_DMAi_pairs = final_results %>%
  distinct(chr, start, end, strand, .keep_all = TRUE) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  group_by(gene_name, pA_type) %>%  # Group by both gene_name and UTR_type
  slice_max(order_by = avg_abs_change_in_usage, n = 1) %>%  # Keep only the row with the highest avg_abs_change_in_usage for each UTR_type within each gene_name
  ungroup() %>%
  group_by(gene_name) %>%
  filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA')) %>%  # Ensure each gene group has at least one 'proximal pA' and one 'distal pA'
  mutate(sign_proximal = sign(max_change_in_usage[pA_type == 'proximal pA']),  # Determine the sign of avg_change_in_usage for proximal
         sign_distal = sign(max_change_in_usage[pA_type == 'distal pA'])) %>%  # Determine the sign of avg_change_in_usage for distal
  filter(sign_proximal != sign_distal) %>%  # Keep only those gene groups where the signs are different
  dplyr::select(-sign_proximal, -sign_distal) %>%  # Optionally remove the temporary sign columns
  ungroup()

#deal with different basemeans for different sites of same gene (because some sites may have been identified in different number of lines)
common_DMAi_pairs = common_DMAi_pairs %>%
  group_by(gene_name) %>%
  mutate(avg_basemean = mean(avg_basemean, na.rm = TRUE)) %>%
  ungroup()

### find regulatable sites not affected by methylation = non-DMAi sites: a) > 0.1 change and padj < 0.05 in any non-DMAi condition
filter_non_DMAi_data <- function(non_DMAi_data_frames, columns_of_interest) {
  # Combine all non-DMAi data frames into one, adding a 'df_index' for identification
  combined_df <- bind_rows(non_DMAi_data_frames, .id = "df_index")
  
  # Filter out rows with NA in change_in_usage or twostep_feature_padj
  filtered_df <- combined_df %>%
    filter(!is.na(change_in_usage) & !is.na(twostep_feature_padj)) %>%
    filter(abs(change_in_usage) > 0.1 & twostep_feature_padj < 0.05)
  
  # Group by columns of interest and check for any valid records, calculate max abs change_in_usage
  results_df <- filtered_df %>%
    group_by(across(all_of(columns_of_interest))) %>%
    summarize(
      max_abs_change_in_usage_non_dmai = max(abs(change_in_usage), na.rm = TRUE),  # Calculate the max absolute change_in_usage
      max_change_in_usage_non_dmai = change_in_usage[which.max(abs(change_in_usage))],  # Get the change_in_usage with the maximum absolute value
      basemean_at_max_change_in_usage = baseMean[which.max(abs(change_in_usage))], # Retrieve basemean from the row with max change_in_usage
      ctrl_usage_at_max_change_in_usage = pA_usage_control[which.max(abs(change_in_usage))], # Retrieve ctrl usage from the row with max change_in_usage
      valid = any(abs(change_in_usage) > 0.1 & twostep_feature_padj < 0.05),
      gene_name = dplyr::first(gene_name),  # Assuming gene_name, feature_id, pA_type do not change within groups
      feature_id = dplyr::first(feature_id),
      pA_type = dplyr::first(pA_type),
      .groups = 'drop'
    ) %>%
    filter(valid)  # Keep only groups where the condition is met
  
  return(results_df)
}

# Apply the function
non_DMAi_results <- filter_non_DMAi_data(non_DMAi_data_frames_with_basemean, columns_of_interest)

# Function to filter these non-DMAi sites to remove sites that also change in DMAi conditions
filter_based_on_DMAi <- function(filtered_non_dmai_results, DMAi_data_frames, columns_of_interest) {
  combined_dmai_df <- bind_rows(DMAi_data_frames, .id = "df_index") %>%
    filter(!is.na(change_in_usage))
  
  # Calculate max abs change_in_usage and average abs change_in_usage per group in DMAi data frames and count occurrences
  dmai_stats <- combined_dmai_df %>%
    group_by(across(all_of(columns_of_interest))) %>%
    summarize(
      max_abs_change_in_usage_dmai = max(abs(change_in_usage), na.rm = TRUE),
      avg_abs_change_in_usage_dmai = mean(abs(change_in_usage), na.rm = TRUE),
      occurrence_count = n_distinct(df_index),
      .groups = 'drop'
    ) %>%
    filter(occurrence_count >= 2)
  
  # Join with filtered_non_dmai_results and apply the filtering criteria
  final_filtered_results <- filtered_non_dmai_results %>%
    inner_join(dmai_stats, by = columns_of_interest) %>%
    filter(
      max_abs_change_in_usage_dmai < 0.67 * max_abs_change_in_usage_non_dmai,
      avg_abs_change_in_usage_dmai < 0.33 * max_abs_change_in_usage_non_dmai,
      avg_abs_change_in_usage_dmai < 0.125
    )
  
  return(final_filtered_results)
}

#run function
filtered_non_DMAi_sites <- filter_based_on_DMAi(non_DMAi_results, DMAi_data_frames, columns_of_interest)

# remove duplicates and retain only pA site pairs (i.e genes that have two sites passing the filters)
non_DMAi_pairs <- filtered_non_DMAi_sites %>%
  # Remove duplicates based on 'chr', 'start', 'end', 'strand'
  distinct(chr, start, end, strand, .keep_all = TRUE) %>%
  # Remove duplicates based on 'feature_id'
  distinct(feature_id, .keep_all = TRUE) %>%
  # Group by 'gene_name'
  group_by(gene_name) %>%
  # Keep only gene_name groups with at least one positive and one negative 'max_change_in_usage_non_dmai'
  filter(any(max_change_in_usage_non_dmai > 0) & any(max_change_in_usage_non_dmai < 0)) %>%
  # For each gene_name, select the row with the highest positive and the row with the lowest negative 'max_change_in_usage_non_dmai'
  do({
    df <- .
    pos_row <- df %>%
      filter(max_change_in_usage_non_dmai > 0) %>%
      slice_max(order_by = max_change_in_usage_non_dmai, n = 1)
    neg_row <- df %>%
      filter(max_change_in_usage_non_dmai < 0) %>%
      slice_min(order_by = max_change_in_usage_non_dmai, n = 1)
    bind_rows(pos_row, neg_row)
  }) %>%
  ungroup() %>%
  # Assign 'pA_type' based on 'strand' and 'start'
  group_by(gene_name) %>%
  mutate(
    strand_gene = unique(strand)[1],
    start_max = max(start),
    pA_type = case_when(
      strand_gene == "+" & start == start_max ~ "distal pA",
      strand_gene == "+" & start != start_max ~ "proximal pA",
      strand_gene == "-" & start == start_max ~ "proximal pA",
      strand_gene == "-" & start != start_max ~ "distal pA"
    )
  ) %>%
  dplyr::select(-strand_gene, -start_max) %>%
  ungroup()

#deal with differing basemeans
non_DMAi_pairs = non_DMAi_pairs %>%
  group_by(gene_name) %>%
  mutate(basemean_at_max_change_in_usage = mean(basemean_at_max_change_in_usage, na.rm = TRUE)) %>%
  ungroup()

##find control sites i.e that are a) on average between -0.075 and 0.075 change in usage, b) don't exceed 0.15 change in usage in any condition, c) have an average usage > 0.15
# Define the columns of interest and thresholds
filter_combined_data <- function(DMAi_list, non_DMAi_list, columns_of_interest) {
  # Step 1: Extract unique combinations from DMAi data frames
  get_unique_sites <- function(data_frames_list, columns_of_interest) {
    unique_sites <- lapply(data_frames_list, function(df) {
      df %>% dplyr::select(all_of(columns_of_interest)) %>% distinct()
    })
    unique_sites_df <- bind_rows(unique_sites) %>% distinct()
    return(unique_sites_df)
  }
  
  DMAi_sites <- get_unique_sites(DMAi_list, columns_of_interest)
  
  # Step 2: Extract unique combinations from non-DMAi data frames
  non_DMAi_sites <- get_unique_sites(non_DMAi_list, columns_of_interest)
  
  # Step 3: Find common combinations between DMAi and non-DMAi sites
  common_sites <- inner_join(DMAi_sites, non_DMAi_sites, by = columns_of_interest)
  
  # Step 4: Merge all data frames from both lists
  merged_list <- c(DMAi_list, non_DMAi_list)
  
  # Combine all data frames into one, adding a 'df_index' for identification
  combined_df <- bind_rows(merged_list, .id = "df_index")
  
  # Step 5: Retain only rows with columns_of_interest present in common_sites
  combined_df <- combined_df %>% semi_join(common_sites, by = columns_of_interest)
  
  # Filter out rows with NA in change_in_usage or mean_usage
  filtered_df <- combined_df %>%
    filter(!is.na(change_in_usage) & !is.na(mean_usage))
  
  # Group by columns of interest and compute necessary statistics
  results_df <- filtered_df %>%
    group_by(across(all_of(columns_of_interest))) %>%
    summarize(
      total_count = n_distinct(df_index),
      max_abs_change_in_usage = max(abs(change_in_usage), na.rm = TRUE),
      avg_abs_change_in_usage = mean(abs(change_in_usage), na.rm = TRUE),
      max_change_in_usage = change_in_usage[which.max(abs(change_in_usage))],
      avg_mean_usage = mean(mean_usage, na.rm = TRUE),
      avg_basemean = mean(baseMean, na.rm = TRUE),
      avg_ctrl_usage = mean(pA_usage_control, na.rm = TRUE),
      gene_name = dplyr::first(gene_name),
      feature_id = dplyr::first(feature_id),
      pA_type = dplyr::first(pA_type),
      .groups = 'drop'
    ) %>%
    # Apply the specified filters
    filter(
      total_count > 1,
      max_abs_change_in_usage < 0.15,
      avg_abs_change_in_usage < 0.075,
      avg_mean_usage > 0.15
    ) %>%
    ungroup()
  
  return(results_df)
}

# Apply the function using your DMAi and non-DMAi lists
common_control_sites = filter_combined_data(
  DMAi_data_frames_with_basemean,
  non_DMAi_data_frames_with_basemean,
  columns_of_interest
)

#as with the DMAi and non-DMAi groups, only keep genes that pairs of pA sites passing the filter
common_control_pairs = common_control_sites %>%  distinct(chr, start, end, strand, .keep_all = TRUE) %>%  
  distinct(chr, start, end, strand, .keep_all = TRUE) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  group_by(gene_name, pA_type) %>%  # Group by both gene_name and UTR_type
  slice_max(order_by = avg_mean_usage, n = 1) %>%  # Keep only the row with the highest avg_abs_change_in_usage for each UTR_type within each gene_name
  ungroup() %>%
  group_by(gene_name) %>%
  filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA')) %>%  # Ensure each gene group has at least one 'proximal pA' and one 'distal pA'
  ungroup()

#deal with different basemeans for different sites of same gene (because some sites may have been identified in different number of lines)
common_control_pairs = common_control_pairs %>%
  group_by(gene_name) %>%
  mutate(avg_basemean = mean(avg_basemean, na.rm = TRUE)) %>%
  ungroup()


### add UTR type info to sites so that we can separate into TUTR and non-TUTR
#read in significant and control csv files containing UTR_type annotation
#read in TUTR and ALE info
# Read each file into a list of dataframes
directory <- "/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/panel_of_cancer_lines/APA_classified/"
DMAi_file_list <- list.files(directory, pattern = ".*sig.*all.*", full.names = TRUE)
DMAi_UTR_type_table_list <- lapply(DMAi_file_list, read.csv)

directory <- "/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig4/non_DMAi_datasets/APA_classified/"
non_DMAi_file_list <- list.files(directory, pattern = ".*sig.*all.*", full.names = TRUE)
non_DMAi_UTR_type_table_list <- lapply(non_DMAi_file_list, read.csv)

directory <- "/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/panel_of_cancer_lines/APA_classified/"
DMAi_ctrl_file_list <- list.files(directory, pattern = ".*ctrl.*all.*", full.names = TRUE)
DMAi_ctrl_UTR_type_table_list <- lapply(DMAi_ctrl_file_list, read.csv)

directory <- "/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig4/non_DMAi_datasets/APA_classified/"
non_DMAi_ctrl_file_list <- list.files(directory, pattern = ".*ctrl.*all.*", full.names = TRUE)
non_DMAi_ctrl_UTR_type_table_list <- lapply(non_DMAi_ctrl_file_list, read.csv)

#combine DMAi and non-DMAi lists
sig_files_list = c(DMAi_UTR_type_table_list,non_DMAi_UTR_type_table_list)
ctrl_files_list = c(DMAi_ctrl_UTR_type_table_list,non_DMAi_ctrl_UTR_type_table_list)

#keep relevant columns
sig_files_list  = lapply(sig_files_list,function(df) {
  df = df[,c('chr','start','end','strand','gene_name',
             'feature_id','twostep_feature_padj','change_in_usage','UTR_type')]
  return(df)
})

ctrl_files_list  = lapply(ctrl_files_list,function(df) {
  df = df[,c('chr','start','end','strand','gene_name',
             'feature_id','twostep_feature_padj','change_in_usage','UTR_type')]
  return(df)
})

#combine all the tables into one table
sig_UTR_type_table = do.call(rbind, sig_files_list)
ctrl_UTR_type_table = do.call(rbind, ctrl_files_list)

#join consistent sites with UTR type tables and then only keep one row per pA site
sig_common_pairs_list = list(common_DMAi_pairs,non_DMAi_pairs)

annotated_sig_sites_list <- lapply(sig_common_pairs_list, function(df) {
  df %>%
    left_join(sig_UTR_type_table %>% dplyr::select(feature_id, UTR_type), by = "feature_id") %>%
    group_by(chr, start, end, strand) %>%
    dplyr::slice(1) %>%
    ungroup()
})

annotated_common_control_pairs  = common_control_pairs  %>%
  left_join(ctrl_UTR_type_table  %>% dplyr::select(feature_id, UTR_type), by = "feature_id") %>%
  group_by(chr, start, end, strand) %>%
  dplyr::slice(1) %>%
  ungroup()

#bring annotated dfs into one list
sites_list = c(annotated_sig_sites_list,list(annotated_common_control_pairs))

#remove pairs with contrasting annotations
sites_list = lapply(sites_list, function(df) {
    df <- df %>%
      group_by(gene_name) %>%
      filter(n_distinct(UTR_type) == 1) %>%
      ungroup()
  return(df)
})

#### isolate TUTR and keep pairs
TUTR_sites_list = lapply(sites_list,function(df){
  df = df %>% filter(UTR_type == 'TUTR')
})

TUTR_sites_list = lapply(TUTR_sites_list,function(df){
  df = df %>% group_by(gene_name) %>%
    filter(n() >= 2)  %>% ungroup()
})

#isolate distal and proximal
TUTR_distal_sites_list = lapply(TUTR_sites_list,function(df){
  df = df %>% filter(pA_type == 'distal pA')
})

TUTR_proximal_sites_list = lapply(TUTR_sites_list,function(df){
  df = df %>% filter(pA_type == 'proximal pA')
})

#find directionality
up_TUTR_distal_sites_list = lapply(TUTR_distal_sites_list,function(df){
  # Check if 'avg_change_in_usage' column is present
  if ("change_in_usage" %in% names(df)) {
    df <- df %>% filter(change_in_usage > 0.075)
  }
  # Else, check if 'max_change_in_usage' column is present
  else if ("max_abs_change_in_usage" %in% names(df)) {
    df <- df %>% filter(max_change_in_usage > 0.075)
  }
  # Additionally check if 'max_change_in_usage_non_DMAi' column is present
  else if ("max_abs_change_in_usage_non_dmai" %in% names(df)) {
    df <- df %>% filter(max_change_in_usage_non_dmai > 0.075)
  }
  return(df)
})

down_TUTR_distal_sites_list = lapply(TUTR_distal_sites_list,function(df){
  # Check if 'change_in_usage' column is present
  if ("change_in_usage" %in% names(df)) {
    df <- df %>% filter(change_in_usage < -0.075)
  }
  # Else, check if 'max_change_in_usage' column is present
  else if ("max_abs_change_in_usage" %in% names(df)) {
    df <- df %>% filter(max_change_in_usage < -0.075)
  }
  # Additionally check if 'max_change_in_usage_non_DMAi' column is present
  else if ("max_abs_change_in_usage_non_dmai" %in% names(df)) {
    df <- df %>% filter(max_change_in_usage_non_dmai < -0.075)
  }
  return(df)
})

up_TUTR_proximal_sites_list = lapply(TUTR_proximal_sites_list,function(df){
  # Check if 'change_in_usage' column is present
  if ("change_in_usage" %in% names(df)) {
    df <- df %>% filter(change_in_usage > 0.075)
  }
  # Else, check if 'max_change_in_usage' column is present
  else if ("max_abs_change_in_usage" %in% names(df)) {
    df <- df %>% filter(max_change_in_usage > 0.075)
  }
  # Additionally check if 'max_change_in_usage_non_DMAi' column is present
  else if ("max_abs_change_in_usage_non_dmai" %in% names(df)) {
    df <- df %>% filter(max_change_in_usage_non_dmai > 0.075)
  }
  return(df)
})

down_TUTR_proximal_sites_list = lapply(TUTR_proximal_sites_list,function(df){
  # Check if 'change_in_usage' column is present
  if ("change_in_usage" %in% names(df)) {
    df <- df %>% filter(change_in_usage < -0.075)
  }
  # Else, check if 'max_change_in_usage' column is present
  else if ("max_abs_change_in_usage" %in% names(df)) {
    df <- df %>% filter(max_change_in_usage < -0.075)
  }
  # Additionally check if 'max_change_in_usage_non_DMAi' column is present
  else if ("max_abs_change_in_usage_non_dmai" %in% names(df)) {
    df <- df %>% filter(max_change_in_usage_non_dmai < -0.075)
  }
  return(df)
})

#make list of common DMAi, non-DMAi, ctrl TUTR distal and proximal sites
all_TUTR_list = c(TUTR_distal_sites_list,TUTR_proximal_sites_list,up_TUTR_distal_sites_list[1:2],
                  down_TUTR_distal_sites_list[1:2],up_TUTR_proximal_sites_list[1:2],down_TUTR_proximal_sites_list[1:2])

### Define function to calculate PA regions
calculate_PA_regions <- function(df) {
  df <- df %>%
    mutate(
      region_start = if_else(strand == "+", end - 100, start - 100),
      region_end   = if_else(strand == "+", end + 100, start + 100)
    )
  return(df)
}

# Define a function to prepare BED columns from a data frame
prepare_bed_df <- function(df) {
  bed_df <- data.frame(
    chr = df$chr,
    start = df$region_start , 
    end = df$region_end,
    gene_name = df$gene_name,
    feature_id = df$feature_id,
    strand = df$strand  
  )
  # Ensure 'start' is at least 0
  bed_df$start <- pmax(bed_df$start, 0)
  return(bed_df)
}

# Apply to each merged data frame
all_TUTR_list_with_pA_regions <- lapply(all_TUTR_list, calculate_PA_regions)

# Apply the function to each data frame in all_TUTR_list_with_pA_regions
all_TUTR_bed_files <- lapply(all_TUTR_list_with_pA_regions, prepare_bed_df)

#write bed files function
write_bed_files <- function(df_list) {
  # Create directories if they don't already exist
  if (!dir.exists("common_DMAi")) {
    dir.create("common_DMAi")
  }
  if (!dir.exists("non_DMAi")) {
    dir.create("non_DMAi")
  }
  
  # Define the mapping of indices to filenames and directories
  mapping <- list(
    # Common DMAi files
    list(index = 1, filename = "common_DMAi_sig_all_TUTR_distal_sites.bed", dir = "common_DMAi"),
    list(index = 4, filename = "common_DMAi_sig_all_TUTR_proximal_sites.bed", dir = "common_DMAi"),
    list(index = 7, filename = "common_DMAi_sig_up_TUTR_distal_sites.bed", dir = "common_DMAi"),
    list(index = 9, filename = "common_DMAi_sig_down_TUTR_distal_sites.bed", dir = "common_DMAi"),
    list(index = 11, filename = "common_DMAi_sig_up_TUTR_proximal_sites.bed", dir = "common_DMAi"),
    list(index = 13, filename = "common_DMAi_sig_down_TUTR_proximal_sites.bed", dir = "common_DMAi"),
    list(index = 3, filename = "control_TUTR_distal_sites.bed", dir = "common_DMAi"),
    list(index = 6, filename = "control_TUTR_proximal_sites.bed", dir = "common_DMAi"),
    # Non-DMAi files
    list(index = 2, filename = "non_DMAi_sig_all_TUTR_distal_sites.bed", dir = "non_DMAi"),
    list(index = 5, filename = "non_DMAi_sig_all_TUTR_proximal_sites.bed", dir = "non_DMAi"),
    list(index = 8, filename = "non_DMAi_sig_up_TUTR_distal_sites.bed", dir = "non_DMAi"),
    list(index = 10, filename = "non_DMAi_sig_down_TUTR_distal_sites.bed", dir = "non_DMAi"),
    list(index = 12, filename = "non_DMAi_sig_up_TUTR_proximal_sites.bed", dir = "non_DMAi"),
    list(index = 14, filename = "non_DMAi_sig_down_TUTR_proximal_sites.bed", dir = "non_DMAi")
  )
  
  # Iterate over the mapping and write each data frame
  for (item in mapping) {
    idx <- item$index
    filename <- item$filename
    dir_name <- item$dir
    full_path <- file.path(dir_name, filename)
    
    # Write the data frame to a BED file
    # Assuming each data frame has the necessary columns in the correct order for a BED file
    write.table(df_list[[idx]], file = full_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

# call function
write_bed_files(all_TUTR_bed_files)

##get gene region bed files
#only want to work with DMAi-lengthened, non-DMAi lengthened and control genes that are > 10kb long
prox_gene_groups = c(down_TUTR_proximal_sites_list[1:2],TUTR_proximal_sites_list[3])
distal_gene_groups = c(up_TUTR_distal_sites_list[1:2],TUTR_distal_sites_list[3])
gene_groups <- Map(rbind, prox_gene_groups, distal_gene_groups)

#read in gencode file (set path to the location of your gencode annotation file)
gencode_gtf_path='gencode.v45.annotation.gtf'

# Import the GTF file
gtf <- import(gencode_gtf_path)

# Filter for genes
genes <- gtf[gtf$type == "gene"]

# Extract metadata columns
genes_df <- as.data.frame(genes)

# Compute TSS based on strand
genes_df$TSS <- ifelse(genes_df$strand == "+", genes_df$start, genes_df$end)

# Define the function to merge a single data frame with genes_df
merge_gene_group <- function(df) {
  # Merge with genes_df on 'gene_name'
  merged_genes <- inner_join(df, genes_df, by = "gene_name")
  # Return the merged data frame
  return(merged_genes)
}

# Apply the merge function to each data frame in gene_group_list
merged_gene_groups <- lapply(gene_groups, merge_gene_group)

#add distance between TSS and prox pA
add_distance_column <- function(df) {
  df %>%
    group_by(gene_name) %>%
    mutate(
      # Extract the start.x coordinate from the row where pA_type == 'proximal pA'
      prox_start = start.x[pA_type == 'proximal pA'][1],
      
      # Calculate TSS_to_pPA_distance using the prox_start value
      TSS_to_pPA_distance = if_else(
        strand.x == '+',
        prox_start - TSS,
        TSS - prox_start
      )
    ) %>%
    ungroup()
}

# Apply the function to each data frame in the list
merged_gene_groups_with_gene_distance <- lapply(merged_gene_groups, add_distance_column)

######make region bed files (for 10kb genes only)

#filter to only keep genes > 10kb
merged_gene_groups_with_gene_distance_tenK = lapply(merged_gene_groups_with_gene_distance,function(df){
  df = df %>% filter(TSS_to_pPA_distance > 10000)
  return(df)
})

###first get 3' UTR region
merged_gene_groups_with_UTR_region <- lapply(
  merged_gene_groups_with_gene_distance_tenK,
  function(df) {
    df %>%
      group_by(gene_name) %>%
      mutate(
        UTR_region_start = if_else(
          strand.x == "+",
          min(end.x) - 1000,   # For plus-strand
          min(start.x) - 1000  # For minus-strand
        ),
        UTR_region_end = if_else(
          strand.x == "+",
          max(end.x) + 1000,   # For plus-strand
          max(start.x) + 1000  # For minus-strand
        )
      ) %>%
      ungroup()
  }
)

#only need one per gene
merged_gene_groups_with_UTR_region_single = lapply(merged_gene_groups_with_UTR_region,function(df){
  df = df %>%
    group_by(gene_name) %>%
    dplyr::slice(1) %>% ungroup()
})

# Define a function to prepare BED columns from a data frame
prepare_bed_df <- function(df) {
  bed_df <- data.frame(
    chr = df$seqnames,
    start = df$UTR_region_start , 
    end = df$UTR_region_end,
    name = df$gene_name,
    score = ".",
    strand = df$strand.x)  
  return(bed_df)
}

# Apply the function to each data frame in gene_groups_with_UTR_regions
UTR_bed_dfs <- lapply(merged_gene_groups_with_UTR_region_single, prepare_bed_df)

#### Define function to calculate upstream of gene regions

#only proximal row per gene
merged_gene_groups_with_gene_distance_tenK_proximal  = lapply(merged_gene_groups_with_gene_distance_tenK,function(df) {
  df = df %>% group_by(gene_name) %>% dplyr::filter(pA_type == 'proximal pA')
  return(df)
})

calculate_upstream_of_gene_regions <- function(df) {
  # Compute regions around TSS based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$TSS - 5000,   # For "+" strand, start is TSS - 5000 
    df$TSS + 2000             # For "-" strand, start is TSS + 2000
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$TSS - 2000,            # For "+" strand, end is TSS - 2000
    df$TSS + 5000             # For "-" strand, end is TSS + 5000
  )
  
  return(df)
}

# Define a function to prepare BED columns from a data frame
prepare_bed_df <- function(df) {
  bed_df <- data.frame(
    chr = df$seqnames,
    start = df$region_start , 
    end = df$region_end,
    name = df$gene_name,
    score = ".",
    strand = df$strand.x  
  )
  # Ensure 'start' is at least 0
  bed_df$start <- pmax(bed_df$start, 0)
  return(bed_df)
}

# Apply to each merged data frame
gene_groups_with_upstream_of_gene_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_proximal, calculate_upstream_of_gene_regions)

# Apply the function to each data frame in gene_groups_with_upstream_of_gene_regions
upstream_of_gene_bed_dfs <- lapply(gene_groups_with_upstream_of_gene_regions, prepare_bed_df)

#### Define function to calculate pre-TSS regions
calculate_pre_tss_regions <- function(df) {
  # Compute regions around TSS based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$TSS - 2000,   # For "+" strand, start is TSS - 2000 
    df$TSS + 1000             # For "-" strand, start is TSS + 1000
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$TSS - 1000,            # For "+" strand, end is TSS - 1000
    df$TSS + 2000             # For "-" strand, end is TSS + 2000
  )
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_pre_tss_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_proximal, calculate_pre_tss_regions)

# Apply the function to each data frame in gene_groups_with_tss_regions
pre_tss_bed_dfs <- lapply(gene_groups_with_pre_tss_regions, prepare_bed_df)

###### Define function to calculate TSS regions
calculate_tss_regions <- function(df) {
  # Calculate regions around TSS
  df$region_start <- pmax(df$TSS - 1000, 1)
  df$region_end <- df$TSS + 1000
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_tss_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_proximal, calculate_tss_regions)

# Apply the function to each data frame in gene_groups_with_tss_regions
tss_bed_dfs <- lapply(gene_groups_with_tss_regions, prepare_bed_df)

### Define function to calculate gene body regions
calculate_gene_body_regions <- function(df) {
  # Compute regions around TSS based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$TSS + 1000,  
    df$start.x + 1000          
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$end.x - 1000,      
    df$TSS - 1000         
  )
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_gene_body_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_proximal, calculate_gene_body_regions)

# Apply the function to each data frame in gene_groups_with_tss_regions
gene_body_bed_dfs <- lapply(gene_groups_with_gene_body_regions, prepare_bed_df)

##get distal pA sites for calculating downstream of gene regions
merged_gene_groups_with_gene_distance_tenK_distal  = lapply(merged_gene_groups_with_gene_distance_tenK,function(df) {
  df = df %>% group_by(gene_name) %>% dplyr::filter(pA_type == 'distal pA')
  return(df)
})

##get post distal pA site regions
calculate_post_dPA_regions <- function(df) {
  # Compute regions around pA based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$end.x + 1000,   
    df$start.x - 2000            
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$end.x + 2000,          
    df$start.x - 1000            
  )
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_post_dPA_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_distal, calculate_post_dPA_regions)

# Apply the function to each data frame in gene_groups_regions
post_dPA_bed_dfs <- lapply(gene_groups_with_post_dPA_regions, prepare_bed_df)

### Define function to calculate downstream of gene region
calculate_downstream_of_gene_regions <- function(df) {
  # Compute regions around pA based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$end.x + 2000,   
    df$start.x - 5000            
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$end.x + 5000,          
    df$start.x - 2000            
  )
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_downstream_of_gene_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_distal, calculate_downstream_of_gene_regions)

# Apply the function to each data frame in gene_groups_regions
downstream_of_gene_bed_dfs <- lapply(gene_groups_with_downstream_of_gene_regions, prepare_bed_df)

#combine all gene region bed files into a single list
gene_group_bed_files = c(upstream_of_gene_bed_dfs,pre_tss_bed_dfs,tss_bed_dfs,gene_body_bed_dfs,UTR_bed_dfs,post_dPA_bed_dfs,downstream_of_gene_bed_dfs)

#write bed
prefixes <- c('DMAi', 'non_DMAi', 'control')
suffixes <- c('_upstream_of_gene.bed','_pre_TSS.bed', '_TSS.bed', '_gene_body.bed',
              '_UTR.bed', '_post_dPA.bed','_downstream_of_gene.bed')

# Loop over the data frames in bed_dfs
for (i in seq_along(gene_group_bed_files)) {
  # Get the current data frame
  bed_df <- gene_group_bed_files[[i]]
  
  # Determine the prefix index (1, 2, or 3)
  prefix_index <- ((i - 1) %% 3) + 1
  prefix <- prefixes[prefix_index]
  
  # Determine the suffix based on group number
  group_number <- ceiling(i / 3)
  suffix <- suffixes[group_number]
  
  # Construct the filename
  filename <- paste0(prefix, suffix)
  
  # Write the data frame to the file
  write.table(bed_df, file = filename, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#### GC boxplot
#Define the directory containing the FASTA files - these need to be generated with getFasta from the bedtools suite, using the bed files generated above as input
#setwd
fasta_dir <- "boxplot_regions/"

# Get the list of FASTA files
fasta_files <- list.files(path = fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

# Initialize an empty list to collect data frames
gc_content_list <- list()

# Define a mapping from file prefixes to condition labels
condition_map <- list(
  'control' = 'control',
  'DMAi' = 'DMAi lengthened',
  'non_DMAi' = 'non-DMAi lengthened'
)

# Loop over each FASTA file
for (fasta_file in fasta_files) {
  # Extract the base name of the file (without path)
  fasta_filename <- basename(fasta_file)
  
  # Remove the .fasta extension
  fasta_filename_no_ext <- sub("\\.fasta$", "", fasta_filename)
  
  # Extract the condition and sequence type using regular expressions
  if (grepl("_upstream_of_gene$", fasta_filename_no_ext)) {
    sequence_type <- 'upstream_of_gene'
    condition_key <- sub("_upstream_of_gene$", "", fasta_filename_no_ext)
  } else if (grepl("_downstream_of_gene$", fasta_filename_no_ext)) {
    sequence_type <- 'downstream_of_gene'
    condition_key <- sub("_downstream_of_gene$", "", fasta_filename_no_ext)
  } else if (grepl("_TSS_to_pPA$", fasta_filename_no_ext)) {
    sequence_type <- 'TSS_to_pPA'
    condition_key <- sub("_TSS_to_pPA$", "", fasta_filename_no_ext)
  } else if (grepl("_pPA_to_dPA$", fasta_filename_no_ext)) {
    sequence_type <- 'pPA_to_dPA'
    condition_key <- sub("_pPA_to_dPA$", "", fasta_filename_no_ext)
  } else if (grepl("_pre_TSS$", fasta_filename_no_ext)) {
    sequence_type <- 'pre_TSS'
    condition_key <- sub("_pre_TSS$", "", fasta_filename_no_ext)
  } else if (grepl("_TSS$", fasta_filename_no_ext)) {
    sequence_type <- 'TSS'
    condition_key <- sub("_TSS$", "", fasta_filename_no_ext)
  } else if (grepl("_gene_body$", fasta_filename_no_ext)) {
    sequence_type <- 'gene_body'
    condition_key <- sub("_gene_body$", "", fasta_filename_no_ext)
  } else if (grepl("_UTR$", fasta_filename_no_ext)) {
    sequence_type <- 'UTR'
    condition_key <- sub("_UTR$", "", fasta_filename_no_ext)
  } else if (grepl("_post_dPA$", fasta_filename_no_ext)) {
    sequence_type <- 'post_dPA'
    condition_key <- sub("_post_dPA$", "", fasta_filename_no_ext)
  } else {
    # Skip files that don't match expected patterns
    next
  }
  
  # Map the condition key to the condition label
  condition_label <- condition_map[[condition_key]]
  if (is.null(condition_label)) {
    condition_label <- condition_key  # Use original key if not in the map
  }
  
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file)
  
  # Calculate GC content for each sequence
  gc_counts <- letterFrequency(sequences, letters = c("G", "C"))
  sequence_lengths <- width(sequences)
  gc_content <- (rowSums(gc_counts) / sequence_lengths) * 100  # GC content as a proportion
  
  # Create a data frame
  df <- data.frame(
    sequence_name = names(sequences),
    gc_content = gc_content,
    condition = condition_label,
    sequence_type = sequence_type,
    stringsAsFactors = FALSE
  )
  
  # Append to the list
  gc_content_list[[length(gc_content_list) + 1]] <- df
}

# Combine all data frames in gc_content_list into one if needed
gc_content_df <- do.call(rbind, gc_content_list)

# Convert 'sequence_type' to factors for controlled plotting
gc_content_df$sequence_type <- factor(gc_content_df$sequence_type, levels = c('upstream_of_gene','pre_TSS','TSS','gene_body','UTR','post_dPA','downstream_of_gene'))

# Create boxplots with the new labels and specified colors
ggplot(gc_content_df, aes(x = sequence_type, y = gc_content, fill = condition)) +
  geom_boxplot() +
  ylab('GC Content (%)') +
  xlab('Gene region') +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 0.85,vjust = 0.9),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),  # Remove legend title
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  scale_fill_manual(
    values = c(
      "DMAi lengthened" = "#e53b3b",
      "non-DMAi lengthened" = "darkgreen",
      "control" = "darkgrey"
    )
  ) +
  scale_x_discrete(labels = c("upstream_of_gene" = "-5kb to -2kb","pre_TSS" = "pre-TSS (-2kb to -1kb)", "TSS" = "TSS",'gene_body' = 'Gene Body','UTR' = "3' UTR",'post_dPA' = "post 3'UTR (+1kb to +2kb)",'downstream_of_gene' = '+2kb to +5kb')) +
  guides(fill = guide_legend(title = NULL))



### create expression matched transcript groups
#make list of lengthened DMAi, non-DMAi and ctrl groups
TUTR_distal_sites_list_for_expression_matching = c(up_TUTR_distal_sites_list[1:2],TUTR_distal_sites_list[3])

# Calculate the percentile threshold of Basemean values in TUTR_sites_list
lower_threshold <- quantile(TUTR_distal_sites_list_for_expression_matching[[1]]$avg_basemean, 0.05, na.rm = TRUE)
upper_threshold_2 <- quantile(TUTR_distal_sites_list_for_expression_matching[[2]]$basemean_at_max_change_in_usage, 0.975, na.rm = TRUE)
upper_threshold_3 <- quantile(TUTR_distal_sites_list_for_expression_matching[[3]]$avg_basemean, 0.8, na.rm = TRUE)

# Filter out rows where Basemean does not meet specified threshold
TUTR_distal_sites_list_for_expression_matching[[1]] <- TUTR_distal_sites_list_for_expression_matching[[1]] %>%
  filter(avg_basemean >= lower_threshold | is.na(avg_basemean))
TUTR_distal_sites_list_for_expression_matching[[2]] <- TUTR_distal_sites_list_for_expression_matching[[2]] %>%
  filter(basemean_at_max_change_in_usage <= upper_threshold_2 | is.na(basemean_at_max_change_in_usage))
TUTR_distal_sites_list_for_expression_matching[[3]] <- TUTR_distal_sites_list_for_expression_matching[[3]] %>%
  filter(avg_basemean <= upper_threshold_3 | is.na(avg_basemean))

#filter proximal tables to only keep genes that are present in newly filtered distal tables
down_TUTR_proximal_sites_list_filtered <- list()

for (i in 1:2) {
  prox_df <- down_TUTR_proximal_sites_list[[i]]
  dist_df <- TUTR_distal_sites_list_for_expression_matching[[i]]
  
  # Filter rows in prox_df where gene_name is present in dist_df
  down_TUTR_proximal_sites_list_filtered[[i]] <- prox_df[prox_df$gene_name %in% dist_df$gene_name, , drop = FALSE]
}

#do same for control table
control_TUTR_proximal_sites_filtered = TUTR_proximal_sites_list[[3]][TUTR_proximal_sites_list[[3]]$gene_name %in% TUTR_distal_sites_list_for_expression_matching[[3]]$gene_name, , drop = FALSE]

#make list of common DMAi, non-DMAi, ctrl TUTR distal and proximal sites
all_expression_matched_TUTR_list = c(TUTR_distal_sites_list_for_expression_matching,down_TUTR_proximal_sites_list_filtered,list(control_TUTR_proximal_sites_filtered))

### Define function to calculate pPA regions
calculate_PA_regions = function(df) {
  df$region_start = df$start - 1000
  df$region_end = df$start + 1000      
  return(df)
}

# Apply to each merged data frame
all_expression_matched_TUTR_with_pA_regions_list <- lapply(all_expression_matched_TUTR_list, calculate_PA_regions)

# Apply the previously made bed file generation function to each data frame in all_expression_matched_TUTR_with_pA_regions_list
all_expression_matched_TUTR_bed_files <- lapply(all_expression_matched_TUTR_with_pA_regions_list, prepare_bed_df)

# write expression matched bed files if required
write_bed_files <- function(df_list) {
  # Create directories if they don't already exist
  if (!dir.exists("common_DMAi")) {
    dir.create("common_DMAi")
  }
  if (!dir.exists("non_DMAi")) {
    dir.create("non_DMAi")
  }
  
  # Define the mapping of indices to filenames and directories
  mapping <- list(
    list(index = 1, filename = "common_DMAi_sig_up_TUTR_distal_sites.bed", dir = "common_DMAi"),
    list(index = 2, filename = "non_DMAi_sig_up_TUTR_distal_sites.bed", dir = "non_DMAi"),
    list(index = 3, filename = "control_TUTR_distal_sites.bed", dir = "common_DMAi"),
    list(index = 4, filename = "common_DMAi_sig_down_TUTR_proximal_sites.bed", dir = "common_DMAi"),
    list(index = 5, filename = "non_DMAi_sig_down_TUTR_proximal_sites.bed", dir = "non_DMAi"),
    list(index = 6, filename = "control_TUTR_proximal_sites.bed", dir = "common_DMAi")
  )
  
  # Iterate over the mapping and write each data frame
  for (item in mapping) {
    idx <- item$index
    filename <- item$filename
    dir_name <- item$dir
    full_path <- file.path(dir_name, filename)
    
    # Write the data frame to a BED file
    # Assuming each data frame has the necessary columns in the correct order for a BED file
    write.table(df_list[[idx]], file = full_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

# call function
write_bed_files(all_expression_matched_TUTR_bed_files)

#### assess baseline proximal pA usage in expression matched groups
DMAi_down_proximal_sites = down_TUTR_proximal_sites_list_filtered[[1]] %>%
  mutate(category = 'DMAi') %>%
  dplyr::rename(prox_usage = avg_ctrl_usage) %>% dplyr::select(prox_usage,category)

non_DMAi_down_proximal_sites = down_TUTR_proximal_sites_list_filtered[[2]] %>%
  mutate(category = 'non-DMAi') %>%
  dplyr::rename(prox_usage = ctrl_usage_at_max_change_in_usage) %>% dplyr::select(prox_usage,category)

ctrl_proximal_sites = control_TUTR_proximal_sites_filtered %>%
  mutate(category = 'control') %>%
  dplyr::rename(prox_usage = avg_ctrl_usage) %>% dplyr::select(prox_usage,category)

#make one combined table and plot
all_sites_table = rbind(DMAi_down_proximal_sites,non_DMAi_down_proximal_sites,ctrl_proximal_sites)

ggplot(all_sites_table, aes(x = category, y = prox_usage, fill = category)) +
  geom_boxplot() +
  ylab("DMSO pPA usage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  theme( 
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 16),                         # Increase axis text (ticks) size
    legend.title = element_text(size = 16),                      # Increase legend title size
    legend.text = element_text(size = 18)                        # Increase legend text size
  )+
  scale_fill_manual(values = c(
    "control" = "darkgrey",  # blue
    "DMAi" = "#e53b3b",  # orange
    "non-DMAi" = "#669b5f"  # green
  ))


# Checking normality for each group
shapiro.test(all_sites_table$prox_usage[all_sites_table$category == "control"])
shapiro.test(all_sites_table$prox_usage[all_sites_table$category == "DMAi"])
shapiro.test(all_sites_table$prox_usage[all_sites_table$category == "non-DMAi"])

ggplot(all_sites_table, aes(x = prox_usage)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(all_sites_table$prox_usage, na.rm = TRUE), 
                                         sd = sd(all_sites_table$prox_usage, na.rm = TRUE)), 
                color = "red", size = 1) +
  facet_wrap(~ category, scales = "free") +  # Separate plots for each category
  labs(title = "Histogram with Normal Curve",
       x = "DMSO pPA usage",
       y = "Density") +
  theme_minimal()

# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(prox_usage ~ category, data = all_sites_table)

# Display the result
kruskal_test_result

# Perform Dunn's test for multiple comparisons
dunn_test_result <- dunn.test(all_sites_table$prox_usage, all_sites_table$category, method = "bonferroni")  #ctrl vs dmai-lengthened: padj =  6.488197e-14, ctrl vs non-dmai-lengthened: padj =3.553605e-24, dmai vs non-dmai-lengthened: padj = 1.251070e-01
