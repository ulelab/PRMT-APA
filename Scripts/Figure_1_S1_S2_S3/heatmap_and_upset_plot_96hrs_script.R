library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
library(UpSetR)

###setwd
setwd("../../Data/Figure_1_S1_S2_S3/APA_CSV_files/96hrs/positional_info/")

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

#restrict to most changed sites for heatmap and only want distal pA rows
most_regulated_sig_sites_list = lapply(sig_sites_list,function(df){
  df = df %>%
    filter(abs(prox_to_dist_shift) > 0.3) %>% filter(pA_type == 'distal pA')
  return(df)
})

#add the sig sites into the original list
most_regulated_sites_list = c(most_regulated_sig_sites_list,all_sites_csv_with_prox_to_dist_shift[4:6])

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
sig_ADMA_most_reg_sites=most_regulated_sites_list[[1]]
sig_SDMA_most_reg_sites=most_regulated_sites_list[[3]]
sig_SDMA_ADMA_most_reg_sites=most_regulated_sites_list[[2]]

all_ADMA_sites=most_regulated_sites_list[[4]]
all_SDMA_sites=most_regulated_sites_list[[6]]
all_SDMA_ADMA_sites=most_regulated_sites_list[[5]]

#merge tables
SDMA_ADMA_most_reg_sites_with_ADMA_data = merge(sig_SDMA_ADMA_most_reg_sites, all_ADMA_sites, by = c('gene_name', 'feature_id', 'chr','start', 'end', 'strand'), all.x = TRUE)
SDMA_ADMA_most_reg_sites_with_ADMA_and_SDMA_data = merge(SDMA_ADMA_most_reg_sites_with_ADMA_data, all_SDMA_sites, by = c('gene_name', 'feature_id', 'chr', 'start', 'end', 'strand'), all.x = TRUE)
# Specify the column names you want to keep
columns_to_keep <- c("gene_name", "feature_id", "chr", "start", "end", 
                     "strand", "UTR_type","ADMAi_change_in_usage", "SDMAi_change_in_usage", "SDMAi_ADMAi_change_in_usage")
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
non_TUTR_sig_sites = all_sig_sites %>% filter(UTR_type != 'TUTR') %>% filter(UTR_type != 'iAPA')

#isolate just distal change in usage for heatmap and multiply by 100 to turn in % change in dPA usage
TUTR_heatmap <- TUTR_sig_sites %>%
  select(8:10) %>%
  mutate(across(everything(), ~ .x * 100))
non_TUTR_heatmap <- non_TUTR_sig_sites %>%
  select(8:10) %>%
  mutate(across(everything(), ~ .x * 100))

#define colour scheme
RdBu=colorRampPalette(rev(brewer.pal(9,'RdBu')))(40)  #extrapolate 9 colours palette out to 45

# Define symmetrical breaks around 0
breaks = seq(-40, 40, length.out = 41)

#make heatmap
pheatmap(non_TUTR_heatmap,clustering_distance_rows = "manhattan",
         clustering_distance_cols = "canberra",color=RdBu,breaks=breaks,show_rownames = FALSE,fontsize_col=14,angle_col='45',legend_breaks = c(-50, -25, 0, 25, 50),   # where numeric labels appear
         legend_labels = c("-50%", "-25%", "0%", "25%", "50%"),cutree_rows = 3,labels_col=c("ADMAi","SDMAi","SDMAi + ADMAi"))

####Upset plot
# ------------------------------------------------------------------------------
# 1) Set working directory and read in files
# ------------------------------------------------------------------------------
setwd("../../Data/Figure_1_S1_S2_S3/APA_CSV_files/96hrs/")

# List all CSV files that have 'sig' in their name
csv_files <- list.files(pattern = "sig")

# Read each CSV file into a list
sig_sites_list <- lapply(csv_files, read_csv)

# Here we assume:
#   sig_sites_list[[1]] == ADMAi
#   sig_sites_list[[2]] == SDMAi_ADMAi
#   sig_sites_list[[3]] == SDMAi

# ------------------------------------------------------------------------------
# 2) Add a 'direction' column to each data frame
# ------------------------------------------------------------------------------
sig_sites_list_with_direction <- lapply(sig_sites_list, function(df) {
  df %>%
    mutate(direction = if_else(
      change_in_usage > 0, 
      "increased", 
      "decreased"
    ))
})

# ------------------------------------------------------------------------------
# 3) Combine feature_id & direction into a single column 'feature_direction'
# ------------------------------------------------------------------------------
sig_sites_list_with_direction <- lapply(sig_sites_list_with_direction, function(df) {
  df %>%
    mutate(feature_direction = paste(feature_id, direction, sep = "_"))
})

# ------------------------------------------------------------------------------
# 4) Name each data frame for clarity & add binary flags
# ------------------------------------------------------------------------------
df_ADMAi <- sig_sites_list_with_direction[[1]] %>%
  mutate(in_ADMAi = TRUE)

df_SDMAi_ADMAi <- sig_sites_list_with_direction[[2]] %>%
  mutate(in_SDMAi_ADMAi = TRUE)

df_SDMAi <- sig_sites_list_with_direction[[3]] %>%
  mutate(in_SDMAi = TRUE)

# ------------------------------------------------------------------------------
# 5) Merge data to identify common 'feature_direction' rows among pairs
#    (Only needed if you want to see which events are in all three.)
# ------------------------------------------------------------------------------
# Merge ADMAi & SDMAi
merge_ADMAi_SDMAi <- merge(df_ADMAi, df_SDMAi, by = "feature_direction")

# Merge ADMAi & SDMAi+ADMAi
merge_ADMAi_SDMAi_ADMAi <- merge(df_ADMAi, df_SDMAi_ADMAi, by = "feature_direction")

# Merge SDMAi & SDMAi+ADMAi
merge_SDMAi_SDMAi_ADMAi <- merge(df_SDMAi, df_SDMAi_ADMAi, by = "feature_direction")

# ------------------------------------------------------------------------------
# 6) Identify genes in all three (common intersection) & label them
# ------------------------------------------------------------------------------
common_all <- merge(merge_ADMAi_SDMAi, merge_ADMAi_SDMAi_ADMAi, by = "feature_direction")
common_all <- merge(common_all, merge_SDMAi_SDMAi_ADMAi, by = "feature_direction")
common_all$in_all_three <- TRUE

# ------------------------------------------------------------------------------
# 7) Create a single master table (all_data) that contains all
#    'feature_direction' plus the presence/absence flags
# ------------------------------------------------------------------------------
all_data <- Reduce(
  f = function(x, y) merge(x, y, by = "feature_direction", all = TRUE),
  x = list(
    df_ADMAi[,       c("feature_direction", "in_ADMAi")],
    df_SDMAi_ADMAi[, c("feature_direction", "in_SDMAi_ADMAi")],
    df_SDMAi[,       c("feature_direction", "in_SDMAi")],
    common_all[,     c("feature_direction", "in_all_three")]
  )
)

# Replace NA with FALSE to indicate absence
all_data[is.na(all_data)] <- FALSE

# ------------------------------------------------------------------------------
# 8) Prepare the data frame for UpSetR
#    - Convert TRUE/FALSE to 1/0
#    - Select columns in the order: ADMAi, SDMAi+ADMAi, SDMAi
# ------------------------------------------------------------------------------
upset_data <- all_data %>%
  mutate(ADMAi         = as.integer(in_ADMAi),
         SDMAi_ADMAi   = as.integer(in_SDMAi_ADMAi),
         SDMAi         = as.integer(in_SDMAi)) %>%
  select(ADMAi, SDMAi_ADMAi, SDMAi)

# (Optional) Rename columns if you want slightly nicer labels
# but it's simplest if these match the sets argument exactly
names(upset_data) <- c("ADMAi", "SDMAi+ADMAi", "SDMAi")

# ------------------------------------------------------------------------------
# 9) Generate the UpSet plot
#    sets = c("ADMAi","SDMAi+ADMAi","SDMAi") must match the column names
# ------------------------------------------------------------------------------
upset(
  upset_data,
  sets = c("SDMAi", "ADMAi", "SDMAi+ADMAi"),
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = c("#F3ADAD","#ED6262", "#9C2A2A"),
  main.bar.color = "black",
  sets.x.label = "APA events",
  text.scale = c(3, 3, 3, 3, 2, 3)
)
