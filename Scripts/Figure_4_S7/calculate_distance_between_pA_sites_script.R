library(tidyverse)
library(data.table)
library(gridExtra)
library(rstatix)

#read in gene group bed files
setwd('../../Data/Figure_4_S7/bed/APA_3seq/200nt_window')

no_mit_distal_sites = fread('no_mit_distal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

no_mit_proximal_sites = fread('no_mit_proximal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

part_mit_distal_sites = fread('DMAi_part_mit_distal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

part_mit_proximal_sites = fread('DMAi_part_mit_proximal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

full_mit_distal_sites = fread('DMAi_full_mit_distal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

full_mit_proximal_sites = fread('DMAi_full_mit_proximal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

#combine tables of same category sites
no_mit_genes = rbind(no_mit_distal_sites,no_mit_proximal_sites)
part_mit_genes = rbind(part_mit_distal_sites,part_mit_proximal_sites)
full_mit_genes = rbind(full_mit_distal_sites,full_mit_proximal_sites)

#make list of genes
all_genes = list(no_mit_genes,part_mit_genes,full_mit_genes)

#make columns calculating pA distance
all_genes_distance_list = lapply(all_genes,function(df){
  df = df %>%
    group_by(gene_name) %>%
    mutate(pA_distance = abs(max(start) - min(start)))
  return(df)
})

#add category column 
all_genes_distance_list[[1]]$category = 'no mitigation'
all_genes_distance_list[[2]]$category = 'DMAi partial mitigation'
all_genes_distance_list[[3]]$category = 'DMAi full mitigation'

#ønly need one row per gene
all_genes_distance_list = lapply(all_genes_distance_list,function(df){
  df = df %>% group_by(gene_name) %>% dplyr::slice(1) %>% ungroup()
  return(df)
})

#combine all the tables into one table
all_genes_with_distance_table = do.call(rbind, all_genes_distance_list)

#refactor category variable
all_genes_with_distance_table$category = factor(all_genes_with_distance_table$category, levels = c("no mitigation","DMAi partial mitigation","DMAi full mitigation"))

#boxplot of 3' UTR length
ggplot(all_genes_with_distance_table, aes(x = category, y = pA_distance, fill = category)) +
  geom_boxplot() +
  ylab("pPA-dPA nt distance") +
  scale_fill_manual(values = c("no mitigation" = "#469de1","DMAi partial mitigation" = "#ff9e4e","DMAi full mitigation" = "#ff594e"))+
  theme_classic() +
  ylim(0,7000) + 
  theme(
    axis.title = element_text(size = 20),                        # Increase axis title size
    axis.text = element_text(size = 17),  
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.75), # Adjusted for readability
    legend.title = element_text(size = 16),                      # Increase legend title size
    legend.text = element_text(size = 14)                        # Increase legend text size
  ) 

# ensure consistent ordering for interpretation
df2 <- all_genes_with_distance_table %>%
  filter(category %in% c("no mitigation", "DMAi full mitigation")) %>%
  mutate(category = droplevels(category))

# Wilcoxon rank-sum (two-sided)
wilcox_res <- df2 %>%
  wilcox_test(pA_distance ~ category) %>%
  adjust_pvalue(method = "BH") %>%   
  add_significance()
wilcox_res
