library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
library(UpSetR)

###setwd
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/3_end_seq_timecourse_LU99/Nobby_APA_analysis/common_atlas/dedup/github/CSV_files/positional_info")

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
sig_sites_list = all_sites_csv_with_prox_to_dist_shift[1:5]

#restrict to most changed sites for heatmap and keep only one row per gene 
most_regulated_sig_sites_list = lapply(sig_sites_list,function(df){
  df = df %>%
    filter(abs(prox_to_dist_shift) > 0.4) %>% filter(pA_type == 'distal pA')
  return(df)
})

#add the sig sites into the original list
most_regulated_sites_list = c(most_regulated_sig_sites_list,all_sites_csv_with_prox_to_dist_shift[6:10])

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
sig_CDK11i_most_reg_sites=most_regulated_sites_list[[1]]
sig_24hrs_most_reg_sites=most_regulated_sites_list[[2]]
sig_48hrs_most_reg_sites=most_regulated_sites_list[[3]]
sig_72hrs_most_reg_sites=most_regulated_sites_list[[4]]
sig_96hrs_most_reg_sites=most_regulated_sites_list[[5]]

all_CDK11i_sites=most_regulated_sites_list[[6]]
all_24hrs_sites=most_regulated_sites_list[[7]]
all_48hrs_sites=most_regulated_sites_list[[8]]
all_72hrs_sites=most_regulated_sites_list[[9]]
all_96hrs_sites=most_regulated_sites_list[[10]]

#merge tables
most_reg_96hrs_sites_with_72_data = merge(sig_96hrs_most_reg_sites, all_72hrs_sites, by = c('gene_name', 'feature_id','chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_96hrs_sites_with_72_and_48_data = merge(most_reg_96hrs_sites_with_72_data, all_48hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_96hrs_sites_with_72_and_48_and_24_data = merge(most_reg_96hrs_sites_with_72_and_48_data, all_24hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_96hrs_sites_with_72_and_48_and_24_and_CDK11_data = merge(most_reg_96hrs_sites_with_72_and_48_and_24_data, all_CDK11i_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
# Specify the column names you want to keep
columns_to_keep <- c("gene_name", "feature_id", "chr", "start", "end", 
                     "strand", "UTR_type","DMAi_24hrs_change_in_usage", "DMAi_48hrs_change_in_usage", "DMAi_72hrs_change_in_usage","DMAi_96hrs_change_in_usage","CDK11i_24hrs_change_in_usage")
most_reg_96hrs_sites_with_72_and_48_and_24_and_CDK11_data = most_reg_96hrs_sites_with_72_and_48_and_24_and_CDK11_data[,columns_to_keep] 

most_reg_72hrs_sites_with_96_data = merge(sig_72hrs_most_reg_sites, all_96hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start' , 'end', 'strand'), all.x = TRUE)
most_reg_72hrs_sites_with_96_and_48_data = merge(most_reg_72hrs_sites_with_96_data, all_48hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_72hrs_sites_with_96_and_48_and_24_data = merge(most_reg_72hrs_sites_with_96_and_48_data, all_24hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_72hrs_sites_with_96_and_48_and_24_and_CDK11_data = merge(most_reg_72hrs_sites_with_96_and_48_and_24_data, all_CDK11i_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_72hrs_sites_with_96_and_48_and_24_and_CDK11_data = most_reg_72hrs_sites_with_96_and_48_and_24_and_CDK11_data[,columns_to_keep]

most_reg_48hrs_sites_with_96_data = merge(sig_48hrs_most_reg_sites, all_96hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_48hrs_sites_with_96_and_72_data = merge(most_reg_48hrs_sites_with_96_data, all_72hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_48hrs_sites_with_96_and_72_and_24_data = merge(most_reg_48hrs_sites_with_96_and_72_data, all_24hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_48hrs_sites_with_96_and_72_and_24_and_CDK11_data = merge(most_reg_48hrs_sites_with_96_and_72_and_24_data, all_CDK11i_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_48hrs_sites_with_96_and_72_and_24_and_CDK11_data = most_reg_48hrs_sites_with_96_and_72_and_24_and_CDK11_data[,columns_to_keep]

most_reg_24hrs_sites_with_96_data = merge(sig_24hrs_most_reg_sites, all_96hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_24hrs_sites_with_96_and_72_data = merge(most_reg_24hrs_sites_with_96_data, all_72hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_24hrs_sites_with_96_and_72_and_48_data = merge(most_reg_24hrs_sites_with_96_and_72_data , all_48hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_24hrs_sites_with_96_and_72_and_48_and_CDK11_data = merge(most_reg_24hrs_sites_with_96_and_72_and_48_data , all_CDK11i_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_24hrs_sites_with_96_and_72_and_48_and_CDK11_data = most_reg_24hrs_sites_with_96_and_72_and_48_and_CDK11_data[,columns_to_keep]

most_reg_CDK11i_sites_with_96_data = merge(sig_CDK11i_most_reg_sites, all_96hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_CDK11i_sites_with_96_and_72_data = merge(most_reg_CDK11i_sites_with_96_data, all_72hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_CDK11i_sites_with_96_and_72_and_48_data = merge(most_reg_CDK11i_sites_with_96_and_72_data , all_48hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_CDK11i_sites_with_96_and_72_and_48_and_24_data = merge(most_reg_CDK11i_sites_with_96_and_72_and_48_data , all_24hrs_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
most_reg_CDK11i_sites_with_96_and_72_and_48_and_24_data = most_reg_CDK11i_sites_with_96_and_72_and_48_and_24_data[,columns_to_keep]

all_sig_sites = rbind(most_reg_96hrs_sites_with_72_and_48_and_24_and_CDK11_data,most_reg_72hrs_sites_with_96_and_48_and_24_and_CDK11_data,most_reg_48hrs_sites_with_96_and_72_and_24_and_CDK11_data,most_reg_24hrs_sites_with_96_and_72_and_48_and_CDK11_data,most_reg_CDK11i_sites_with_96_and_72_and_48_and_24_data) %>%
  distinct(feature_id, .keep_all = TRUE) %>% distinct(gene_name, .keep_all = TRUE)

#remove rows where there is an NA mean diff value
all_sig_sites = all_sig_sites %>%
  filter(rowSums(is.na(.)) == 0)

#isolate TUTR and non TUTR
TUTR_sig_sites = all_sig_sites %>% filter(UTR_type == 'TUTR')
non_TUTR_sig_sites = all_sig_sites %>% filter(UTR_type != 'TUTR')

#isolate just distal change in usage for heatmap and multiply by 100 to turn in % change in dPA usage
TUTR_heatmap <- TUTR_sig_sites %>%
  select(8:12) %>%
  mutate(across(everything(), ~ .x * 100))
non_TUTR_heatmap <- non_TUTR_sig_sites %>%
  select(8:12) %>%
  mutate(across(everything(), ~ .x * 100))

#define colour scheme
RdBu=colorRampPalette(rev(brewer.pal(9,'RdBu')))(40)  #extrapolate 9 colours palette out to 45

# Define symmetrical breaks around 0
breaks = seq(-40, 40, length.out = 41)

#make heatmap
pheatmap(non_TUTR_heatmap,clustering_distance_rows = "manhattan",
         clustering_distance_cols = "canberra",color=RdBu,breaks=breaks,show_rownames = FALSE,fontsize_col=14,angle_col='45',legend_breaks = c(-50, -25, 0, 25, 50),   # where numeric labels appear
         legend_labels = c("-50%", "-25%", "0%", "25%", "50%"),cutree_rows = 2,labels_col=c("DMAi 24hrs","DMAi 48hrs","DMAi 72hrs","DMAi 96hrs","CDK11i/TOPKi 24hrs"))
