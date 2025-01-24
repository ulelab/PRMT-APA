library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
library(UpSetR)

###setwd - change to location of directory containing the polyA sites with positional info CSVs (generated from the drimseq.R and classification_of_APA_events.R R scripts)
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig1/ADMAi_vs_SDMAi_vs_double/positional_info")

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
sig_sites_list = all_sites_csv_with_prox_to_dist_shift[1:3]

#restrict to most changed sites for heatmap and keep only one row per gene 
most_regulated_sig_sites_list = lapply(sig_sites_list,function(df){
  df = df %>%
    filter(abs(prox_to_dist_shift) > 0.4)
  return(df)
})

#add the sig sites into the original list
most_regulated_sites_list = c(most_regulated_sig_sites_list,all_sites_csv_with_prox_to_dist_shift[4:6])

#retain only necessary columns and one row per gene
most_regulated_sites_list = lapply(most_regulated_sites_list, function(df) {
  column_names_to_keep <- c("gene_name", "feature_id", "chr", "start",
                            "end", "strand", "pA_type", "UTR_type",
                            "condition", "prox_to_dist_shift")
  
  # Check if the specified columns exist in the dataframe
  existing_columns <- column_names_to_keep[column_names_to_keep %in% names(df)]
  
  # Subset the dataframe to include only existing specified columns
  df <- df[, existing_columns, drop = FALSE] %>% group_by(gene_name) %>%
    dplyr::slice(1) %>% ungroup()
  return(df)
})

#add condition onto prox_to_dist shift column 
most_regulated_sites_list = lapply(most_regulated_sites_list, function(df) {
  # Create a new column name based on the first entry in the 'condition' column
  new_colname = paste0(df$condition[1], "_prox_to_dist_shift")
  
  # Check if 'prox_to_dist_shift' exists in the dataframe
  if ("prox_to_dist_shift" %in% names(df)) {
    # Change the name of the 'prox_to_dist_shift' column
    names(df)[names(df) == "prox_to_dist_shift"] = new_colname
  }
  
  # Return the modified dataframe
  return(df)
})

#extract tables
sig_ADMA_most_reg_sites=most_regulated_sites_list[[1]]
sig_SDMA_most_reg_sites=most_regulated_sites_list[[2]]
sig_SDMA_ADMA_most_reg_sites=most_regulated_sites_list[[3]]

all_ADMA_sites=most_regulated_sites_list[[4]]
all_SDMA_sites=most_regulated_sites_list[[5]]
all_SDMA_ADMA_sites=most_regulated_sites_list[[6]]

#merge tables
SDMA_ADMA_most_reg_sites_with_ADMA_data = merge(sig_SDMA_ADMA_most_reg_sites, all_ADMA_sites, by = c('gene_name', 'feature_id', 'chr','start', 'end', 'strand'), all.x = TRUE)
SDMA_ADMA_most_reg_sites_with_ADMA_and_SDMA_data = merge(SDMA_ADMA_most_reg_sites_with_ADMA_data, all_SDMA_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)

# Specify the column names you want to keep
columns_to_keep <- c("gene_name", "feature_id", "chr", "start", "end", 
                     "strand", "UTR_type","ADMAi_prox_to_dist_shift", "SDMAi_prox_to_dist_shift", "SDMAiADMAi_prox_to_dist_shift")
SDMA_ADMA_most_reg_sites_with_ADMA_and_SDMA_data=SDMA_ADMA_most_reg_sites_with_ADMA_and_SDMA_data[, columns_to_keep]

ADMA_most_reg_sites_with_SDMA_data = merge(sig_ADMA_most_reg_sites, all_SDMA_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
ADMA_most_reg_sites_with_SDMA_and_SDMA_ADMA_data = merge(ADMA_most_reg_sites_with_SDMA_data, all_SDMA_ADMA_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
ADMA_most_reg_sites_with_SDMA_and_SDMA_ADMA_data = ADMA_most_reg_sites_with_SDMA_and_SDMA_ADMA_data[,columns_to_keep]

SDMA_most_reg_sites_with_ADMA_data = merge(sig_SDMA_most_reg_sites, all_ADMA_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
SDMA_most_reg_sites_with_ADMA_and_SDMA_ADMA_data = merge(SDMA_most_reg_sites_with_ADMA_data, all_SDMA_ADMA_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
SDMA_most_reg_sites_with_ADMA_and_SDMA_ADMA_data = SDMA_most_reg_sites_with_ADMA_and_SDMA_ADMA_data[,columns_to_keep]

all_sig_sites = rbind(SDMA_ADMA_most_reg_sites_with_ADMA_and_SDMA_data,ADMA_most_reg_sites_with_SDMA_and_SDMA_ADMA_data,SDMA_most_reg_sites_with_ADMA_and_SDMA_ADMA_data) %>%
  distinct(feature_id, .keep_all = TRUE) %>% distinct(gene_name, .keep_all = TRUE)

#remove rows where there is an NA mean diff value
all_sig_sites = all_sig_sites %>%
  filter(rowSums(is.na(.)) == 0)

#isolate TUTR and non TUTR
TUTR_sig_sites = all_sig_sites %>% filter(UTR_type == 'TUTR')
non_TUTR_sig_sites = all_sig_sites %>% filter(UTR_type != 'TUTR')

#isolate just prox-to-dist shifts for heatmap
TUTR_heatmap = TUTR_sig_sites[,c(8:10)]
non_TUTR_heatmap = non_TUTR_sig_sites[,c(8:10)]

#define colour scheme
RdBu=colorRampPalette(rev(brewer.pal(9,'RdBu')))(30)  #extrapolate 9 colours palette out to 30

# Define symmetrical breaks around 0
breaks = seq(-1, 1, length.out = 31)

#make heatmap
pheatmap(non_TUTR_heatmap,clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",color=RdBu,breaks=breaks,fontsize_row = 2,fontsize_col=14,angle_col='45',cutree_rows = 2,labels_col=c("ADMAi","SDMAi","SDMAi + ADMAi"))

####Upset plot
###setwd
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig1/ADMAi_vs_SDMAi_vs_double/")

#### read in files
# List all CSV files
csv_files <- list.files(pattern = "sig")
# Read each CSV file into a list
sig_sites_list <- lapply(csv_files, read_csv)

#add direction
sig_sites_list_with_direction = lapply(sig_sites_list,function(df){
  df = df %>% mutate(direction = if_else(
    change_in_usage > 0,
    'increased',
    'decreased'
  ))
  return(df)
})

#combine feature_id and direction data for upset
sig_sites_list_with_direction = lapply(sig_sites_list_with_direction,function(df){
  df = df %>%
    mutate(feature_direction = paste(feature_id, direction, sep = "_"))
  return(df)
})

# Merge SDMA and ADMA
merge_SDMAi_and_ADMAi <- merge(sig_sites_list_with_direction[[1]], sig_sites_list_with_direction[[2]], by = "feature_direction", suffixes = c(".SDMAi", ".ADMAi"))

# Merge SDMA and ADMA+SDMA 
merge_SDMAi_and_ADMAi_SDMAi <- merge(sig_sites_list_with_direction[[1]], sig_sites_list_with_direction[[3]], by = "feature_direction", suffixes = c(".SDMAi", ".SDMAi_ADMAi"))

# Merge SDMA and ADMA+SDMA 
merge_ADMAi_and_ADMAi_SDMAi <- merge(sig_sites_list_with_direction[[2]], sig_sites_list_with_direction[[3]], by = "feature_direction", suffixes = c(".ADMAi", ".SDMAi_ADMAi"))


# Identify genes that are in all datasets with specified direction criteria
common_all <- merge(merge_SDMAi_and_ADMAi,
                    merge_SDMAi_and_ADMAi_SDMAi, by = "feature_direction")
common_all <- merge(common_all,
                    merge_ADMAi_and_ADMAi_SDMAi, by = "feature_direction")
common_all$in_all_three <- TRUE

# Create binary flags for presence in the datasets according to the criteria
sig_sites_list_with_direction[[1]]$in_SDMAi <- TRUE
sig_sites_list_with_direction[[2]]$in_ADMAi <- TRUE
sig_sites_list_with_direction[[3]]$in_SDMAi_ADMAi <- TRUE

# Combine all information
all_data <- Reduce(function(x, y) merge(x, y, by = "feature_direction", all = TRUE), 
                   list(sig_sites_list_with_direction[[1]][, c("feature_direction", "in_SDMAi")], 
                        sig_sites_list_with_direction[[2]][, c("feature_direction", "in_ADMAi")], 
                        sig_sites_list_with_direction[[3]][, c("feature_direction", "in_SDMAi_ADMAi")],
                        common_all[, c("feature_direction", "in_all_three")]))
all_data[is.na(all_data)] <- FALSE

# Convert data to appropriate format for UpSetR
upset_data <- all_data %>%
  mutate(SDMAi = as.integer(in_SDMAi),
         ADMAi = as.integer(in_ADMAi),
         SDMAi_ADMAi = as.integer(in_SDMAi_ADMAi)) %>%
  select(SDMAi, ADMAi, SDMAi_ADMAi)

# Renaming columns to more descriptive names
names(upset_data) <- c("SDMAi", "ADMAi", "SDMAi+ADMAi")

# Generate the UpSet plot
upset(upset_data, 
      sets = c("SDMAi", "ADMAi", "SDMAi+ADMAi"),
      order.by = "freq",
      keep.order = TRUE,
      sets.bar.color = c("#F3ADAD","#ED6262","#9C2A2A"),
      main.bar.color = "black",
      sets.x.label = "APA events",
      text.scale = c(3, 3, 3, 3, 2, 3))





