library(tidyverse)
library(data.table)
library(gridExtra)
library(rstatix)

##DMAi vs non-DMAi sites
#read in gene group bed files
setwd('/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/bed/test_triple_removal/1000/')

DMAi_lengthened_distal_sites = fread('common_DMAi/common_DMAi_sig_up_TUTR_distal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

DMAi_lengthened_proximal_sites = fread('common_DMAi/common_DMAi_sig_down_TUTR_proximal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

ctrl_distal_sites = fread('common_DMAi/control_TUTR_distal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

ctrl_proximal_sites = fread('common_DMAi/control_TUTR_proximal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

non_DMAi_lengthened_distal_sites = fread('non_DMAi/non_DMAi_sig_up_TUTR_distal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

non_DMAi_lengthened_proximal_sites = fread('non_DMAi/non_DMAi_sig_down_TUTR_proximal_sites.bed') %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    gene_name = V4,
    feature_ID = V5,
    strand = V6
  )

#combine tables of same category sites
DMAi_lengthened_genes = rbind(DMAi_lengthened_distal_sites,DMAi_lengthened_proximal_sites)
non_DMAi_lengthened_genes = rbind(non_DMAi_lengthened_distal_sites,non_DMAi_lengthened_proximal_sites)
ctrl_genes = rbind(ctrl_distal_sites,ctrl_proximal_sites)

#make list of genes
all_genes = list(DMAi_lengthened_genes,non_DMAi_lengthened_genes,ctrl_genes)

#make columns calculating pA distance - do this by substracting the highest start value - lowest start value (even though these 1000 nt window files, this cancels itself out during the calculation, i.e (distal pA - 500) - (proximal pA - 500) = distal pA - proximal pA 
all_genes_distance_list = lapply(all_genes,function(df){
  df = df %>%
    group_by(gene_name) %>%
    mutate(pA_distance = abs(max(start) - min(start)))
  return(df)
})

#add category column 
all_genes_distance_list[[1]]$category = 'DMAi'
all_genes_distance_list[[2]]$category = 'non-DMAi'
all_genes_distance_list[[3]]$category = 'control'

#ønly need one row per gene
all_genes_distance_list = lapply(all_genes_distance_list,function(df){
  df = df %>% group_by(gene_name) %>% dplyr::slice(1) %>% ungroup()
  return(df)
})

#combine all the tables into one table
all_genes_with_distance_table = do.call(rbind, all_genes_distance_list)

#boxplot of 3' UTR length
ggplot(all_genes_with_distance_table, aes(x = category, y = pA_distance, fill = category)) +
  geom_boxplot() +
  ylab("pPAS-dPAS nt distance") +
  ylim(0, 6500) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  theme( 
    axis.title = element_text(size = 20), 
    axis.text = element_text(size = 16),                         # Increase axis text (ticks) size
    legend.title = element_text(size = 16),                      # Increase legend title size
    legend.text = element_text(size = 20)                        # Increase legend text size
  )+
  scale_fill_manual(values = c(
    "control" = "darkgrey",  
    "DMAi" = "#ff594e",  
    "non-DMAi" = "#1A85FF"  
  ))

# ensure consistent ordering for interpretation
all_genes_with_distance_table <- all_genes_with_distance_table %>%
  mutate(category = factor(category, levels = c("control","DMAi","non-DMAi")))

# pairwise Wilcoxon vs control 
pw_dist <- all_genes_with_distance_table %>%
  pairwise_wilcox_test(pA_distance ~ category,
                       ref.group = "control",
                       p.adjust.method = "BH") %>%  
  filter(group2 %in% c("DMAi","non-DMAi")) %>%
  mutate(comparison = paste(group1, "vs", group2))


