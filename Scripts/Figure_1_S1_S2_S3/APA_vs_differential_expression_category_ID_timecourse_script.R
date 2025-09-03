library(tidyverse)
library(purrr)

#setwd
setwd("/Users/k2362866/Documents/AZ_postdoc/3_end_seq_timecourse_LU99/Nobby_APA_analysis/common_atlas/dedup/github/TUTR_vs_ALE/")

TUTR_24hrs = read.csv('all_DMAi_24hrs_sites_with_UTR_type_info.csv') %>% filter(UTR_type == 'TUTR')
TUTR_48hrs = read.csv('all_DMAi_48hrs_sites_with_UTR_type_info.csv') %>% filter(UTR_type == 'TUTR')
TUTR_72hrs = read.csv('all_DMAi_72hrs_sites_with_UTR_type_info.csv') %>% filter(UTR_type == 'TUTR')
TUTR_96hrs = read.csv('all_DMAi_96hrs_sites_with_UTR_type_info.csv') %>% filter(UTR_type == 'TUTR')

#make list
TUTR_list = list(TUTR_24hrs,TUTR_48hrs,TUTR_72hrs,TUTR_96hrs)

#add prox to dist shift column 
TUTR_list = lapply(TUTR_list,function(df){
  df = df %>% group_by(gene_name) %>%
    mutate(prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] - change_in_usage[pA_type == "proximal pA"])
  return(df)
})


##read in whole gene DESeq2 files
#setwd
setwd("/Users/k2362866/Documents/AZ_postdoc/3_end_seq_timecourse_LU99/DESeq2_github")

DE_table_24hrs = read.csv('DESeq2_results_24hrs_all_genes.csv')
DE_table_48hrs = read.csv('DESeq2_results_48hrs_all_genes.csv')
DE_table_72hrs = read.csv('DESeq2_results_72hrs_all_genes.csv')
DE_table_96hrs = read.csv('DESeq2_results_96hrs_all_genes.csv')

#make list
DE_table_list = list(DE_table_24hrs,DE_table_48hrs,DE_table_72hrs,DE_table_96hrs)

#change colname
DE_table_list = lapply(DE_table_list,function(df){
  colnames(df)[7]='gene_name'
  df = df[,c(3,6,7)]
  return(df)
})

#merge tables in list
# Specify the common column
common_column = "gene_name"

# Merge data frames from list1 and list2 based on the common column
merged_df_list = lapply(seq_along(TUTR_list), function(i) {
  # Extract data frames from list1 and list2
  df1 = TUTR_list[[i]]
  df2 = DE_table_list[[i]]
  
  # Merge based on the common column
  merged_df <- merge(df1, df2, by = common_column, all.x = TRUE)
  
  return(merged_df)
})

#change NAs to 1
merged_df_list_no_NA = lapply(merged_df_list, function(df) {
  df$twostep_feature_padj[is.na(df$twostep_feature_padj)] <- 1
  return(df)
})

#only need one row per gene
merged_df_proximal_list = lapply(merged_df_list_no_NA,function(df){
  df = df %>% filter(pA_type == 'proximal pA')
  return(df)
})

#add APA_and_gene_class column
merged_df_list_with_APA_and_gene_class = lapply(merged_df_proximal_list, function(df) {
  df = df %>% mutate(APA_and_gene_class = case_when(
    log2FoldChange >= 0.585 & padj <= 0.05 & prox_to_dist_shift >= 0.2 ~ "3' UTR lengthening + gene upregulated",
    log2FoldChange >= 0.585 & padj <= 0.05 & prox_to_dist_shift <= 0.2 ~ "no 3' UTR lengthening + gene upregulated",
    log2FoldChange <= -0.585 & padj <= 0.05 & prox_to_dist_shift >= 0.2 ~ "3' UTR lengthening + gene downregulated",
    log2FoldChange <= -0.585 & padj <= 0.05 & prox_to_dist_shift <= 0.2 ~ "no 3' UTR lengthening + gene downregulated",
    log2FoldChange >= -0.585 & log2FoldChange <= 0.585 & prox_to_dist_shift >= 0.2 ~ "3' UTR lengthening + no expression change",
    TRUE ~ "No change"
  ))
  return(df)
})

#remove NAs
merged_df_list_with_APA_and_gene_class <- lapply(merged_df_list_with_APA_and_gene_class, function(df) {
  df <- df %>% filter(!is.na(APA_and_gene_class))
  return(df)
})

#factorise category
# Specify the factor levels
factor_levels <- c("3' UTR lengthening + gene upregulated",
                   "3' UTR lengthening + no expression change",
                   "3' UTR lengthening + gene downregulated",
                   "no 3' UTR lengthening + gene upregulated",
                   "no 3' UTR lengthening + gene downregulated",
                   'No change')

# Apply factor level specification
merged_df_list_with_APA_and_gene_class = lapply(merged_df_list_with_APA_and_gene_class,function(df){
  df$APA_and_gene_class = factor(df$APA_and_gene_class, levels = factor_levels)
  return(df)
})

# Define improved color scheme for categories
custom_colors <- c(
  "3' UTR lengthening + gene upregulated" = "#670303",   
  "3' UTR lengthening + gene downregulated" = "#f98a87",   
  "3' UTR lengthening + no expression change" = "#009E73",    
  "no 3' UTR lengthening + gene upregulated" = "#004278",    
  "no 3' UTR lengthening + gene downregulated" = "#5dbfed",   
  'No change' = "#808080"                                        
)

# Create the plot
ggplot(merged_df_list_with_APA_and_gene_class[[4]], 
       aes(x = log2FoldChange, y = prox_to_dist_shift, 
           color = APA_and_gene_class, size = APA_and_gene_class)) +
  geom_point() +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(0.2), linetype = "dashed", color = "black") +
  labs(
    x = "transcript log2FoldChange",
    y = "proximal to distal usage shift",
    color = "APA and DE Class",  # Set the title for the color legend
    size = "APA and DE Class"    # Set the title for the size legend if applicable
  ) +
  theme_classic() +
  xlim(-2.5,2.5) +
  ylim(-1,1) +
  scale_color_manual(values = custom_colors) +
  scale_size_manual(values = c(rep(1, 5), 0.5)) +  # Smaller size for 'No change'
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

#isolate jernej gene groups
up_and_lengthened = merged_df_list_with_APA_and_gene_class[[4]] %>% filter(APA_and_gene_class == "3' UTR lengthening + gene upregulated")   
down_and_lengthened = merged_df_list_with_APA_and_gene_class[[4]] %>% filter(APA_and_gene_class ==  "3' UTR lengthening + gene downregulated")
lengthened_no_change = merged_df_list_with_APA_and_gene_class[[4]] %>% filter(APA_and_gene_class == "3' UTR lengthening + no expression change")  %>% filter(prox_to_dist_shift > 0.4) #apply further filter to get the gene groups the same size
up_not_lengthened = merged_df_list_with_APA_and_gene_class[[4]] %>% filter(APA_and_gene_class == "no 3' UTR lengthening + gene upregulated")
down_not_lengthened = merged_df_list_with_APA_and_gene_class[[4]] %>% filter(APA_and_gene_class == "no 3' UTR lengthening + gene downregulated")

#write to file
#setwd
setwd("/Users/k2362866/Documents/AZ_postdoc/3_end_seq_timecourse_LU99/Nobby_APA_analysis/common_atlas/dedup/github/APA_vs_DE/gene_groups_for_metascape/96hrs")

write_csv(up_not_lengthened,'up_not_lengthened.csv')
write_csv(down_not_lengthened,'down_not_lengthened.csv')
write_csv(up_and_lengthened,'up_and_lengthened.csv')
write_csv(down_and_lengthened,'down_and_lengthened.csv')
write_csv(lengthened_no_change,'lengthened_no_change_0_4_filter.csv')

#turn into list
metascape_group_list = list(up_not_lengthened,down_not_lengthened,up_and_lengthened,down_and_lengthened,lengthened_no_change)

# 1 ── make sure your list elements are named -------------------------------
names(metascape_group_list) <- c(
  "up_not_lengthened",
  "down_not_lengthened",
  "up_and_lengthened",
  "down_and_lengthened",
  "lengthened_no_change"
)

# 2 ── build the five filtered tables ---------------------------------------
full_tables <- imap(
  metascape_group_list,
  ~ {
    # .x = one element of the list (df or vector); .y = its name
    genes <- if (is.data.frame(.x)) .x$gene_name else .x
    
    merged_df_list[[4]]  %>%
      filter(gene_name %in% genes) %>%      # keep matching rows
      mutate(APA_and_gene_class = .y)
  }
)

#only keep relevant columns 
full_tables <- lapply(full_tables, function(df) {
  df <- df[, c(1:3,13:16,20,25)]
  return(df)
})

#write to file
setwd("/Users/k2362866/Documents/AZ_postdoc/files_for_charlotte/github/metascape_groups/")
write_csv(full_tables[[1]],'up_not_lengthened.csv')
write_csv(full_tables[[2]],'down_not_lengthened.csv')
write_csv(full_tables[[3]],'up_and_lengthened.csv')
write_csv(full_tables[[4]],'down_and_lengthened.csv')
write_csv(full_tables[[5]],'lengthened_no_change.csv')

#change gene class column name
full_tables = lapply(full_tables,function(df){
  colnames(df)[9]='APA_and_gene_class_96hrs'
  return(df)
})

## add DE changes to each table from every timepoint
# order must match the list order in `DE_table_list`
time_points <- c("24hrs", "48hrs", "72hrs", "96hrs")

## create a helper list that already has the right column names
DE_ready <- imap(DE_table_list, ~
                   .x %>%
                   dplyr::select(gene_name, log2FoldChange) %>%
                   dplyr::rename(!!paste0("log2FoldChange_", time_points[.y]) := log2FoldChange)
)

## loop over every full-table, tacking on the four columns
full_tables <- map(
  full_tables,
  ~ purrr::reduce(
    DE_ready,
    dplyr::left_join,
    .init = .x,   # start with the current full_table
    by = "gene_name"
  )
)

## add prox to dist shifts from each timepoint
## 1.  Prepare the four tiny tables that hold just the two columns we need
PDS_ready <- imap(merged_df_proximal_list, ~
                    .x %>%                                   # the i-th DF in the list
                    dplyr::select(gene_name, prox_to_dist_shift) %>%
                    dplyr::rename(
                      !!paste0("prox_to_dist_shift_", time_points[.y]) :=
                        prox_to_dist_shift
                    )
)

## 2.  Add those four columns to every full-table
full_tables <- purrr::map(
  full_tables,
  ~ purrr::reduce(
    PDS_ready,
    dplyr::left_join,
    .init = .x,           # start with the current full_table
    by    = "gene_name"
  )
)

#combine into one big list
big_table = bind_rows(full_tables)

#write to file
setwd("/Users/k2362866/Documents/AZ_postdoc/files_for_charlotte/github/metascape_groups/")
write_csv(big_table,'all_96hrs_GO_groups_with_LFC_and_APA.csv')

## create bed files
# 1. pick only the rows you want from merged_df_list_no_NA[[4]]
rows_to_add <- merged_df_list_no_NA[[4]] %>%          # source table
  semi_join(up_and_lengthened, by = "gene_name") %>%  # keep genes already in up_and_lengthened
  filter(pA_type == "distal pA") %>%                  # …and only distal pA rows
  mutate(APA_and_gene_class = "3' UTR lengthened + gene upregulated") %>%  # add the extra column
  select(names(up_and_lengthened))                    # make sure column order matches

# 2. bind them onto the original table
up_and_lengthened_prox_and_dist <- bind_rows(up_and_lengthened, rows_to_add)

#do same for down_and_lengthened
rows_to_add <- merged_df_list_no_NA[[4]] %>%          # source table
  semi_join(down_and_lengthened, by = "gene_name") %>%  # keep genes already in down_and_lengthened
  filter(pA_type == "distal pA") %>%                  # …and only distal pA rows
  mutate(APA_and_gene_class = "3' UTR lengthened + gene downregulated") %>%  # add the extra column
  select(names(down_and_lengthened))                    # make sure column order matches

# 2. bind them onto the original table
down_and_lengthened_prox_and_dist <- bind_rows(down_and_lengthened, rows_to_add)

# make list
sites_list = list(up_and_lengthened_prox_and_dist,down_and_lengthened_prox_and_dist)

#create windows for STREME etc
#add new columns defining window for fasta file
APA_sites_with_window_list = lapply(sites_list, function(df) {
  df <- df %>%
    mutate(
      region_start = if_else(strand == "+", end - 1000, start - 1000),
      region_end   = if_else(strand == "+", end + 1000, start + 1000)
    )
  return(df)
})

#proximal sites
proximal_sites_with_window_list = lapply(APA_sites_with_window_list,function(df){
  df = df %>% dplyr::filter(pA_type == 'proximal pA')
  return(df)
})

#distal sites
distal_sites_with_window_list = lapply(APA_sites_with_window_list,function(df){
  df = df %>% dplyr::filter(pA_type == 'distal pA')
  return(df)
})

all_APA_sites_with_window_list = c(distal_sites_list,proximal_sites_list)

#keep relevant columns
APA_sites_with_window_list = lapply(APA_sites_with_window_list,function(df){
  df = df %>% select(c("chr", "region_start", "region_end",
                       "gene_name", "feature_id", "strand"))
  return(df)
})

#write bed
setwd("/Users/k2362866/Documents/AZ_postdoc/3_end_seq_timecourse_LU99/Nobby_APA_analysis/common_atlas/dedup/github/APA_vs_DE/bed/96hrs/2000")

write.table(all_APA_sites_with_window_list[[1]], file = 'up_and_lengthened_distal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[2]], file = 'down_and_lengthened_distal_sites.bed',
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[3]], file = 'up_and_lengthened_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[4]], file = 'down_and_lengthened_proximal_sites.bed',
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##create expression matched groups
#read in deseq2 file
whole_gene_96hrs_DE = read.csv("/Users/k2362866/Documents/AZ_postdoc/3_end_seq_timecourse_LU99/DESeq2_github/DESeq2_results_96hrs_all_genes.csv")
colnames(whole_gene_96hrs_DE)[7] = 'gene_name'

#add basemean to gene group tables 
# Loop over each data frame in the list and perform left_join with whole_gene_DE
gene_group_with_basemean_list <- lapply(sites_list, function(df) {
  # Perform the left join to add Basemean column based on matching gene_name
  left_join(df, whole_gene_96hrs_DE %>% dplyr::select(gene_name, baseMean), by = "gene_name")
})

# Calculate the percentile threshold of Basemean values in gene_group_list
lower_threshold <- quantile(gene_group_with_basemean_list[[1]]$baseMean, 0.14, na.rm = TRUE)
upper_threshold_2 <- quantile(gene_group_with_basemean_list[[2]]$baseMean,0.86, na.rm = TRUE)
upper_threshold_3 <- quantile(gene_group_with_basemean_list[[1]]$baseMean,0.99, na.rm = TRUE)

# Filter out rows where Basemean does not meet specified threshold
gene_group_with_basemean_list[[1]] <- gene_group_with_basemean_list[[1]] %>%
  dplyr::filter(baseMean >= lower_threshold | is.na(baseMean))
gene_group_with_basemean_list[[2]] <- gene_group_with_basemean_list[[2]] %>%
  dplyr::filter(baseMean <= upper_threshold_2 | is.na(baseMean))
gene_group_with_basemean_list[[1]] <- gene_group_with_basemean_list[[1]] %>%
  dplyr::filter(baseMean <= upper_threshold_3 | is.na(baseMean))

#add new columns defining window for fasta file
APA_sites_with_window_list = lapply(gene_group_with_basemean_list, function(df) {
  df <- df %>%
    mutate(
      region_start = if_else(strand == "+", end - 1000, start - 1000),
      region_end   = if_else(strand == "+", end + 1000, start + 1000)
    )
  return(df)
})

#isolate distal and proximal
distal_sites_list = lapply(APA_sites_with_window_list,function(df){
  df = df %>% dplyr::filter(pA_type == 'distal pA')
})

proximal_sites_list = lapply(APA_sites_with_window_list,function(df){
  df = df %>% dplyr::filter(pA_type == 'proximal pA')
})

all_APA_sites_with_window_list = c(distal_sites_list,proximal_sites_list)

#keep relevant columns
all_APA_sites_with_window_list = lapply(all_APA_sites_with_window_list,function(df){
  df = df %>% select(c("chr", "region_start", "region_end",
                       "gene_name", "feature_id", "strand"))
  return(df)
})


#setwd and write to file
setwd("/Users/k2362866/Documents/AZ_postdoc/3_end_seq_timecourse_LU99/Nobby_APA_analysis/common_atlas/dedup/github/APA_vs_DE/bed/96hrs/2000/expression_matched/")

#write bed
write.table(all_APA_sites_with_window_list[[1]], file = 'up_and_lengthened_distal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[2]], file = 'down_and_lengthened_distal_sites.bed',
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[3]], file = 'up_and_lengthened_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[4]], file = 'down_and_lengthened_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

### check proximal pA usage for expression matched groups
#only keep proximal pA rows and proximal pA usage + category columns
proximal_pA_usage_data_list = lapply(gene_group_with_basemean_list[1:2],function(df){
  df = df %>% dplyr::filter(pA_type == 'proximal pA') %>% dplyr::select(pA_usage_control,APA_and_gene_class) 
  return(df)
})

#combine into one dataframe for plotting
proximal_pA_usage_data = rbind(proximal_pA_usage_data_list[[1]],proximal_pA_usage_data_list[[2]])

#ggplot
ggplot(proximal_pA_usage_data, aes(x = APA_and_gene_class, y = pA_usage_control, fill = APA_and_gene_class)) +
  geom_boxplot() +
  ylab("baseline proximal pA usage") +
  ylim(0,0.9) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13),                        # Increase axis title size
    axis.text = element_text(size = 12),  
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.75), # Adjusted for readability
    legend.title = element_text(size = 12),                      # Increase legend title size
    legend.text = element_text(size = 11)                        # Increase legend text size
  ) 
