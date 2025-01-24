library(tidyverse)
library(data.table)

##read in significant TUTR site files
#read in sig DMAi cell lines TUTR sites - change directory to the location of your TUTR significant site csv files are (generated from classification_of_APA_event_script.R script)
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/panel_of_cancer_lines/APA_classified")
DMAi_data_frames = list.files(pattern = "^sig.*\\.csv")

# List to store data frames
sig_DMAi_TUTR_csv_files = list()

# Loop through each CSV file and read it
for (file in DMAi_data_frames) {
  df <- read_csv(file)
  sig_DMAi_TUTR_csv_files[[file]] <- df
}

#read in non-DMAi condition TUTR sites
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig3/non_DMAi_datasets/APA_classified")
non_DMAi_data_frames = list.files(pattern = "^sig.*\\.csv")

# List to store data frames
sig_non_DMAi_TUTR_csv_files = list()

# Loop through each CSV file and read it
for (file in non_DMAi_data_frames) {
  df <- read_csv(file)
  sig_non_DMAi_TUTR_csv_files[[file]] <- df
}

#combine lists 
sig_TUTR_dfs = c(sig_DMAi_TUTR_csv_files,sig_non_DMAi_TUTR_csv_files)

## read in ctrl TUTR files
#read in ctrl DMAi cell lines TUTR sites - - change directory to the location of your TUTR ctrl site csv files are (generated from TUTR_vs_ALE_and_jitterplot_script.R script)
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/panel_of_cancer_lines/APA_classified")
DMAi_data_frames = list.files(pattern = "^ctrl.*\\.csv")

# List to store data frames
ctrl_DMAi_TUTR_csv_files = list()

# Loop through each CSV file and read it
for (file in DMAi_data_frames) {
  df <- read_csv(file)
  ctrl_DMAi_TUTR_csv_files[[file]] <- df
}

#read in non-DMAi condition TUTR sites
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/TUTR_vs_ALE/ctrl_CSV/non_DMAi/")
non_DMAi_data_frames = list.files(pattern = "^ctrl.*\\.csv")

# List to store data frames
ctrl_non_DMAi_TUTR_csv_files = list()

# Loop through each CSV file and read it
for (file in non_DMAi_data_frames) {
  df <- read_csv(file)
  ctrl_non_DMAi_TUTR_csv_files[[file]] <- df
}

#combine lists 
ctrl_TUTR_dfs = c(ctrl_DMAi_TUTR_csv_files,ctrl_non_DMAi_TUTR_csv_files)

#for motif analysis, only retain the most changing pA sites (and we want to keep only pairs that pass this threshold)
most_sig_TUTR_dfs = lapply(sig_TUTR_dfs,function(df){
  df = df %>% filter(abs(change_in_usage) > 0.15) %>%
    group_by(gene_name) %>%
    filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA')) %>%  # Ensure each gene group has at least one 'proximal pA' and one 'distal pA'
    mutate(sign_proximal = sign(change_in_usage[pA_type == 'proximal pA']),  # Determine the sign of avg_change_in_usage for proximal
           sign_distal = sign(change_in_usage[pA_type == 'distal pA'])) %>%  # Determine the sign of avg_change_in_usage for distal
    filter(sign_proximal != sign_distal) %>%  # Keep only those gene groups where the signs are different
    select(-sign_proximal, -sign_distal) %>%  # Optionally remove the temporary sign columns
    ungroup() 
})

#for control sites, we want to make sure both sites in pair to be used 
both_used_ctrl_TUTR_dfs = lapply(ctrl_TUTR_dfs,function(df){
  df = df %>% filter(mean_usage > 0.15) %>% 
    group_by(gene_name) %>%
    filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA'))  %>% ungroup()
  return(df)
})

## isolate distal and proximal
sig_TUTR_distal_sites_list = lapply(most_sig_TUTR_dfs,function(df){
  df = df %>% filter(pA_type == 'distal pA')
})

sig_TUTR_proximal_sites_list = lapply(most_sig_TUTR_dfs,function(df){
  df = df %>% filter(pA_type == 'proximal pA')
})

ctrl_TUTR_distal_sites_list = lapply(both_used_ctrl_TUTR_dfs,function(df){
  df = df %>% filter(pA_type == 'distal pA')
})

ctrl_TUTR_proximal_sites_list = lapply(both_used_ctrl_TUTR_dfs,function(df){
  df = df %>% filter(pA_type == 'proximal pA')
})

#find directionality of sig sites
up_TUTR_distal_sites_list = lapply(sig_TUTR_distal_sites_list,function(df){
  df = df %>% filter(change_in_usage > 0)
  return(df)
})

down_TUTR_distal_sites_list = lapply(sig_TUTR_distal_sites_list,function(df){
  df = df %>% filter(change_in_usage < 0)
  return(df)
})

up_TUTR_proximal_sites_list = lapply(sig_TUTR_proximal_sites_list,function(df){
  df = df %>% filter(change_in_usage > 0)
  return(df)
})

down_TUTR_proximal_sites_list = lapply(sig_TUTR_proximal_sites_list,function(df){
  df = df %>% filter(change_in_usage < 0)
  return(df)
})

#make list of common DMAi, non-DMAi, ctrl TUTR distal and proximal sites
all_TUTR_list = c(sig_TUTR_distal_sites_list,sig_TUTR_proximal_sites_list,up_TUTR_distal_sites_list,
                  down_TUTR_distal_sites_list,up_TUTR_proximal_sites_list,down_TUTR_proximal_sites_list,
                  ctrl_TUTR_distal_sites_list,ctrl_TUTR_proximal_sites_list)

### Define function to calculate pPA regions
calculate_PA_regions = function(df) {
  df$region_start = df$start - 100
  df$region_end = df$start + 100          
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

# Apply the function to each data frame in gene_groups_with_tss_regions
all_TUTR_bed_files <- lapply(all_TUTR_list_with_pA_regions, prepare_bed_df)

# Function to write data frames to BED files with specified filenames and directories e.g.
setwd('../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/bed')

# Define the prefixes and suffixes
prefixes <- c('GP2D', 'H838', 'HCT116_p53minus', 'HCT116_p53plus', 'LU99', 'MCF7', 'PANC0403', 'SKMEL2', 'SUM149PT', 'U2OS', 'CBPi', 'CDK11i', 'CDK12i', 'CLKi_C', 'CLKi_N', 'HNRNPC_KD', 'HUB_N_vs_T_DMSO', 'JTE607', 'METTL3i', 'siCFIM25', 'siELAVL1', 'siPCF11', 'siSam68', 'siSUPT4H1', 'siXRN2')

suffixes <- c("_sig_all_TUTR_distal_sites.bed", "_sig_all_TUTR_proximal_sites.bed", "_sig_up_TUTR_distal_sites.bed", "_sig_down_TUTR_distal_sites.bed", "_sig_up_TUTR_proximal_sites.bed", "_sig_down_TUTR_proximal_sites.bed", "_control_TUTR_distal_sites.bed", "_control_TUTR_proximal_sites.bed")

# Loop through each data frame in your list
for (i in 1:length(all_TUTR_bed_files)) {
  
  # Determine the suffix and prefix indices
  suffix_index <- ceiling(i / 25)
  prefix_index <- ((i - 1) %% 25) + 1
  
  # Get the corresponding prefix and suffix
  prefix <- prefixes[prefix_index]
  suffix <- suffixes[suffix_index]
  
  # Access the current data frame
  df <- all_TUTR_bed_files[[i]]
  
  # Create the directory if it doesn't exist
  dir.create(prefix, showWarnings = FALSE)
  
  # Construct the filename and filepath
  filename <- paste0(prefix, suffix)
  filepath <- file.path(prefix, filename)
  
  # Write the data frame to a bed file
  write.table(df, file = filepath, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
