library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)

###setwd
setwd("../../Data/Figure_2_S4/APA_CSV_files/panel_of_cancer_lines/positional_info")

#### read in files
# List all CSV files
csv_files <- list.files(pattern = "\\.csv$")
# Read each CSV file into a list
all_sites_csv <- lapply(csv_files, read_csv)

#add proximal to distal shift column
all_sites_csv_with_prox_to_dist_shift <- lapply(all_sites_csv,function(df) {
  df <- df %>%
    group_by(gene_name) %>%
    mutate(prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] - change_in_usage[pA_type == "proximal pA"]) %>%
    ungroup()
  return(df)
})

#separate into sig sites and all sites for the subsequent filtering step
sig_sites_list = all_sites_csv_with_prox_to_dist_shift[1:10]

#restrict to most changed sites for heatmap and keep only one row per gene 
most_regulated_sig_sites_list = lapply(sig_sites_list,function(df){
  df = df %>%
    dplyr::filter(abs(prox_to_dist_shift) > 0.3) %>% dplyr::filter(pA_type == 'distal pA')
  return(df)
})

#add the sig sites into the original list
most_regulated_sites_list = c(most_regulated_sig_sites_list,all_sites_csv_with_prox_to_dist_shift[11:20])

#retain only necessary columns and one row per gene
most_regulated_sites_list = lapply(most_regulated_sites_list, function(df) {
  column_names_to_keep <- c("gene_name", "feature_id", "chr", "start",
                            "end", "strand", "pA_type", "UTR_type",
                            "condition", "change_in_usage")
  
  # Check if the specified columns exist in the dataframe
  existing_columns <- column_names_to_keep[column_names_to_keep %in% names(df)]
  
  # Subset the dataframe to include only existing specified columns
  df <- df[, existing_columns, drop = FALSE] %>% group_by(gene_name) %>%
    dplyr::slice(1) %>% ungroup()
  return(df)
})

#add condition onto change_in_usage column 
most_regulated_sites_list = lapply(most_regulated_sites_list, function(df) {
  # Create a new column name based on the first entry in the 'condition' column
  new_colname = paste0(df$condition[1], "_change_in_usage")
  
  # Check if 'change_in_usage' exists in the dataframe
  if ("change_in_usage" %in% names(df)) {
    # Change the name of the 'change_in_usage' column
    names(df)[names(df) == "change_in_usage"] = new_colname
  }
  
  # Return the modified dataframe
  return(df)
})


#extract tables
sig_GP2D_sites=most_regulated_sites_list[[1]]
sig_H838_sites=most_regulated_sites_list[[2]]
sig_HCTminus_sites=most_regulated_sites_list[[3]]
sig_HCTplus_sites=most_regulated_sites_list[[4]]
sig_LU99_sites=most_regulated_sites_list[[5]]
sig_MCF7_sites=most_regulated_sites_list[[6]]
sig_PANC0403_sites=most_regulated_sites_list[[7]]
sig_SKMEL2_sites=most_regulated_sites_list[[8]]
sig_SUM149PT_sites=most_regulated_sites_list[[9]]
sig_U2OS_sites=most_regulated_sites_list[[10]]

all_GP2D_sites=most_regulated_sites_list[[11]]
all_H838_sites=most_regulated_sites_list[[12]]
all_HCTminus_sites=most_regulated_sites_list[[13]]
all_HCTplus_sites=most_regulated_sites_list[[14]]
all_LU99_sites=most_regulated_sites_list[[15]]
all_MCF7_sites=most_regulated_sites_list[[16]]
all_PANC0403_sites=most_regulated_sites_list[[17]]
all_SKMEL2_sites=most_regulated_sites_list[[18]]
all_SUM149PT_sites=most_regulated_sites_list[[19]]
all_U2OS_sites=most_regulated_sites_list[[20]]

# Define the list of conditions
conditions <- c("GP2D", "H838", "HCT116_p53minus", "HCT116_p53plus", "LU99",
                "MCF7", "PANC0403", "SKMEL2", "SUM149PT", "U2OS")

# Assuming 'most_regulated_sites_list' contains the dataframes
# The first 10 elements are significant sites dataframes
sig_sites_list <- most_regulated_sites_list[1:10]
# The next 10 elements are all sites dataframes
all_sites_list <- most_regulated_sites_list[11:20]

# Name the lists by condition for easy access
names(sig_sites_list) <- conditions
names(all_sites_list) <- conditions

# Adjusted function to process each condition
process_condition <- function(condition_name, sig_sites_list, all_sites_list, conditions) {
  # Get the significant sites dataframe for the current condition
  sig_df <- sig_sites_list[[condition_name]]
  
  # Ensure 'change_in_usage' column is named correctly
  shift_col_name <- paste0(condition_name, "_change_in_usage")
  if ("change_in_usage" %in% names(sig_df)) {
    sig_df <- sig_df %>% rename_with(~ shift_col_name, .cols = "change_in_usage")
  }
  
  # Select necessary columns from the significant sites dataframe
  sig_df <- sig_df %>%
    dplyr::select(gene_name, feature_id, chr, start, end, strand, UTR_type, all_of(shift_col_name))
  
  # Get the conditions excluding the current one
  other_conditions <- setdiff(conditions, condition_name)
  
  # Initialize the merged dataframe with the significant sites dataframe
  merged_df <- sig_df
  
  # Iterate over other conditions and merge their data
  for (other_cond in other_conditions) {
    # Get the all sites dataframe for the other condition
    all_df <- all_sites_list[[other_cond]]
    
    # Ensure 'change_in_usage' column is named correctly
    other_shift_col_name <- paste0(other_cond, "_change_in_usage")
    if ("change_in_usage" %in% names(all_df)) {
      all_df <- all_df %>% rename_with(~ other_shift_col_name, .cols = "change_in_usage")
    }
    
    # Select necessary columns from the all sites dataframe
    all_df <- all_df %>%
      dplyr::select(gene_name, feature_id, chr, start, end, strand, all_of(other_shift_col_name))
    
    # Merge the current merged_df with the all_df of the other condition
    merged_df <- merge(merged_df, all_df, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
  }
  
  # After merging, select the desired columns
  shift_columns <- grep("change_in_usage", names(merged_df), value = TRUE)
  columns_to_keep <- c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand', 'UTR_type', shift_columns)
  merged_df <- merged_df[, columns_to_keep]
  
  return(merged_df)
}

# Apply the process_condition function to each condition
processed_dfs <- lapply(conditions, function(condition_name) {
  process_condition(condition_name, sig_sites_list, all_sites_list, conditions)
})

# Combine all processed dataframes into one
all_sig_sites <- bind_rows(processed_dfs) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  distinct(gene_name, .keep_all = TRUE)

# Remove rows with NA values in shift columns
# Identify all shift columns dynamically
shift_columns <- grep("change_in_usage", names(all_sig_sites), value = TRUE)
all_sig_sites_complete_data <- all_sig_sites %>%
  filter(rowSums(is.na(select(., all_of(shift_columns)))) == 0)

#isolate TUTR and non TUTR
TUTR_sig_sites = all_sig_sites_complete_data %>% filter(UTR_type == 'TUTR')
non_TUTR_sig_sites = all_sig_sites_complete_data %>% filter(UTR_type != 'TUTR')

#isolate just distal change in usage for heatmap and multiply by 100 to turn in % change in dPA usage
TUTR_heatmap <- TUTR_sig_sites %>%
  select(8:17) %>%
  mutate(across(everything(), ~ .x * 100))
non_TUTR_heatmap <- non_TUTR_sig_sites %>%
  select(8:17) %>%
  mutate(across(everything(), ~ .x * 100))

#define colour scheme
RdBu=colorRampPalette(rev(brewer.pal(9,'RdBu')))(40)  #extrapolate 9 colours palette out to 45

# Define symmetrical breaks around 0
breaks = seq(-35, 35, length.out = 41)

#make heatmap
pheatmap(TUTR_heatmap,clustering_distance_rows = "manhattan",
         clustering_distance_cols = "canberra",color=RdBu,breaks=breaks,show_rownames = FALSE,fontsize_col=14,angle_col='45',legend_breaks = c(-50, -25, 0, 25, 50),   # where numeric labels appear
         legend_labels = c("-50%", "-25%", "0%", "25%", "50%"),cutree_rows = 2,labels_col=c("GP2D","NCI-H838","HCT116 p53 -/-","HCT116 p53 +/+","LU-99","MCF7","PANC0403","SK-MEL-2","SUM149PT","U2OS"))

