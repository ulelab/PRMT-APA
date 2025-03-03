library(ggplot2)
library(tidyverse)

#setwd
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/dedup/github/TUTR_vs_ALE/CSV")

#Read in pA site files and only retain TUTR genes
TUTR_ADMAi = read.csv('all_ADMAi_sites_with_UTR_type_info.csv') %>% filter(UTR_type == 'TUTR')
TUTR_SDMAi = read.csv('all_SDMAi_sites_with_UTR_type_info.csv') %>% filter(UTR_type == 'TUTR')
TUTR_SDMAiADMAi = read.csv('all_SDMAi_ADMAi_sites_with_UTR_type_info.csv') %>% filter(UTR_type == 'TUTR')

#make list
TUTR_list = list(TUTR_ADMAi,TUTR_SDMAi,TUTR_SDMAiADMAi)

#add prox to dist shift column 
TUTR_list = lapply(TUTR_list,function(df){
  df = df %>% group_by(gene_name) %>%
  mutate(prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] - change_in_usage[pA_type == "proximal pA"])
  return(df)
})

##read in whole gene DESeq2 files
#setwd
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/DESeq2_github")

DE_table_ADMAi = read.csv('3_end_seq_DMSO_vs_ADMAi_96hrs_LFC_DE_all_genes.csv')
DE_table_SDMAi = read.csv('3_end_seq_DMSO_vs_SDMAi_96hrs_LFC_DE_all_genes.csv')
DE_table_ADMAiSDMAi = read.csv('3_end_seq_DMSO_vs_SDMAi_ADMAi_96hrs_LFC_DE_all_genes.csv')

#make list
DE_table_list = list(DE_table_ADMAi,DE_table_SDMAi,DE_table_ADMAiSDMAi)

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
merged_df_list = lapply(merged_df_list, function(df) {
  df$twostep_feature_padj[is.na(df$twostep_feature_padj)] <- 1
  return(df)
})

#only need one row per gene
merged_df_list = lapply(merged_df_list,function(df){
  df = df %>% filter(pA_type == 'proximal pA')
  return(df)
})

#pull out SDMAi_ADMAi and define categories
df_SDMAi_ADMAi = merged_df_list[[3]] %>%  mutate(APA_and_gene_class = case_when(
  log2FoldChange > 0.585 & padj < 0.05 & prox_to_dist_shift > 0.2 ~ "3' UTR lengthening + gene upregulated",
  log2FoldChange > 0.585 & padj < 0.05 & (twostep_feature_padj > 0.05 | prox_to_dist_shift < 0.2) ~ "no 3' UTR lengthening + gene upregulated",
  log2FoldChange < -0.585 & padj < 0.05 & prox_to_dist_shift > 0.2 ~ "3' UTR lengthening + gene downregulated",
  log2FoldChange < -0.585 & padj < 0.05 & (twostep_feature_padj > 0.05 | prox_to_dist_shift < 0.2) ~ "no 3' UTR lengthening + gene downregulated",
  log2FoldChange > -0.585 & log2FoldChange < 0.585 & prox_to_dist_shift > 0.2 ~ "3' UTR lengthening + no expression change",
  TRUE ~ "No change"
))

#remove NAs
df_SDMAi_ADMAi = df_SDMAi_ADMAi %>% filter(!is.na(APA_and_gene_class))

#factorise category
# Specify the factor levels
factor_levels <- c("3' UTR lengthening + gene upregulated",
                   "3' UTR lengthening + no expression change",
                   "3' UTR lengthening + gene downregulated",
                   "no 3' UTR lengthening + gene upregulated",
                   "no 3' UTR lengthening + gene downregulated",
                   'No change')

# Apply factor level specification
df_SDMAi_ADMAi$APA_and_gene_class <- factor(df_SDMAi_ADMAi$APA_and_gene_class, levels = factor_levels)

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
ggplot(df_SDMAi_ADMAi, 
       aes(x = log2FoldChange, y = prox_to_dist_shift, 
           color = APA_and_gene_class, size = APA_and_gene_class)) +
  geom_point() +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(0.2), linetype = "dashed", color = "black") +
  labs(
    x = "gene expression log2FoldChange",
    y = "proximal to distal usage shift",
    color = "APA and DE Class",  # Set the title for the color legend
    size = "APA and DE Class"    # Set the title for the size legend if applicable
  ) +
  theme_classic() +
  xlim(-2.25,2.25) +
  ylim(-1,1) +
  scale_color_manual(values = custom_colors) +
  scale_size_manual(values = c(rep(1.2, 5), 0.6)) +  # Smaller size for 'No change'
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )


#isolate gene groups
up_and_lengthened = df_SDMAi_ADMAi %>% filter(APA_and_gene_class == "3' UTR lengthening + gene upregulated")   
down_and_lengthened = df_SDMAi_ADMAi %>% filter(APA_and_gene_class ==  "3' UTR lengthening + gene downregulated")
lengthened_no_change = df_SDMAi_ADMAi %>% filter(APA_and_gene_class == "3' UTR lengthening + no expression change")  %>% filter(prox_to_dist_shift > 0.35 & twostep_gene_padj < 0.05) #apply further filter to get the gene groups the same size
up_not_lengthened = df_SDMAi_ADMAi %>% filter(APA_and_gene_class == "no 3' UTR lengthening + gene upregulated")
down_not_lengthened = df_SDMAi_ADMAi %>% filter(APA_and_gene_class == "no 3' UTR lengthening + gene downregulated")

#write to file
#setwd
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/dedup/github/APA_vs_DE/gene_groups_for_metascape/SDMAi_ADMAi/jernej_method")

write_csv(up_not_lengthened,'up_not_lengthened.csv')
write_csv(down_not_lengthened,'down_not_lengthened.csv')
write_csv(up_and_lengthened,'up_and_lengthened.csv')
write_csv(down_and_lengthened,'down_and_lengthened.csv')
write_csv(lengthened_no_change,'lengthened_no_change_0_35_filter_0_05sig.csv')

