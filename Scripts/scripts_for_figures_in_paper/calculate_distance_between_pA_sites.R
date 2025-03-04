library(tidyverse)
library(data.table)
library(gridExtra)
library(dunn.test)

##DMAi vs non-DMAi sites
#read in gene group bed files
setwd('/Users/llywelyngriffith/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/bed/test_triple_removal/1000/')

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
    "DMAi" = "#e53b3b",  
    "non-DMAi" = "#669b5f"  
  ))

# Checking normality for each group
shapiro.test(all_genes_with_distance_table$pA_distance[all_genes_with_distance_table$category == "control"])
shapiro.test(all_genes_with_distance_table$pA_distance[all_genes_with_distance_table$category == "DMAi"])
shapiro.test(all_genes_with_distance_table$pA_distance[all_genes_with_distance_table$category == "non-DMAi"])

ggplot(all_genes_with_distance_table, aes(x = pA_distance)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(all_genes_with_distance_table$pA_distance, na.rm = TRUE), 
                                         sd = sd(all_genes_with_distance_table$pA_distance, na.rm = TRUE)), 
                color = "red", size = 1) +
  facet_wrap(~ category, scales = "free") +  # Separate plots for each category
  labs(title = "Histogram with Normal Curve",
       x = "pPAS-dPAS nt distance",
       y = "Density") +
  theme_minimal()

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(pA_distance ~ category, data = all_genes_with_distance_table)

# Print the test result
print(kruskal_test_result)

# Perform Dunn's test for multiple comparisons
dunn_test_result <- dunn.test(all_genes_with_distance_table$pA_distance, all_genes_with_distance_table$category, method = "bonferroni") #ctrl vs dmai-lengthened: padj = 8.717070e-27, ctrl vs non-dmai-lengthened: padj = 1.437234e-13, dmai vs non-dmai-lengthened: padj = 1.494149e-04


#####perform same analysis on siCFIM25 category groups
#read in gene group bed files
setwd('/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/200/')

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

# Checking normality for each group
shapiro.test(all_genes_with_distance_table$pA_distance[all_genes_with_distance_table$category == "no mitigation"])
shapiro.test(all_genes_with_distance_table$pA_distance[all_genes_with_distance_table$category == "DMAi partial mitigation"])
shapiro.test(all_genes_with_distance_table$pA_distance[all_genes_with_distance_table$category == "DMAi full mitigation"])

ggplot(all_genes_with_distance_table, aes(x = pA_distance)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(all_genes_with_distance_table$pA_distance, na.rm = TRUE), 
                                         sd = sd(all_genes_with_distance_table$pA_distance, na.rm = TRUE)), 
                color = "red", size = 1) +
  facet_wrap(~ category, scales = "free") +  # Separate plots for each category
  labs(title = "Histogram with Normal Curve",
       x = "pPAS-dPAS nt distance",
       y = "Density") +
  theme_minimal()

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(pA_distance ~ category, data = all_genes_with_distance_table)

# Print the test result
print(kruskal_test_result)

# Perform Dunn's test for multiple comparisons
dunn_test_result <- dunn.test(all_genes_with_distance_table$pA_distance, all_genes_with_distance_table$category, method = "bonferroni") #no mitigation vs full mitigation padj =  

dunn_test_result
