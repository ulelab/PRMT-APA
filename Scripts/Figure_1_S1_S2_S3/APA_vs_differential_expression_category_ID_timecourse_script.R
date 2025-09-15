library(tidyverse)
library(purrr)

#setwd
setwd("../../Data/Figure_1_S1_S2_S3/APA_CSV_files/timecourse
/APA_classified/")

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
setwd("../../Data/Figure_1_S1_S2_S3/DESeq2_tables/timecourse/")

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
setwd("../../Data/Figure_1_S1_S2_S3/metascape_gene_lists/")

write_csv(up_not_lengthened,'up_not_lengthened.csv')
write_csv(down_not_lengthened,'down_not_lengthened.csv')
write_csv(up_and_lengthened,'up_and_lengthened.csv')
write_csv(down_and_lengthened,'down_and_lengthened.csv')
write_csv(lengthened_no_change,'lengthened_no_change_0_4_filter.csv')
