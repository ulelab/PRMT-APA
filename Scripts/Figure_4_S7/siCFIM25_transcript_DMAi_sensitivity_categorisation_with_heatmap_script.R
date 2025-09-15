library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(ggpattern)
library(GenomicFeatures)
library(dunn.test)      

#setwd
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/CSV_files/positional_info/siCFIM25")

# List all CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")  # This regex matches files that end with .csv

# Read each CSV file and store them in a list called 'siCFIM25_files'
siCFIM25_files <- lapply(csv_files, function(file) {
  df <- read.csv(file)  # Or use read_csv(file) for faster loading with readr
  
  # Select necessary columns
  df <- df %>%
    dplyr::select("gene_name", "feature_id", "chr", "start", "end", "strand",
           "pA_usage_control", "pA_usage_treatment", "change_in_usage", "prox_to_dist_shift",
           "twostep_feature_padj", "condition", "pA_type", "mean_usage")
  
  # Assume 'condition' is the same across all rows, take the first one
  condition_value <- unique(df$condition)[1]
  
  # List of columns to rename
  columns_to_rename <- c("prox_to_dist_shift", "twostep_feature_padj", "change_in_usage", 
                         "mean_usage", "pA_usage_control", "pA_usage_treatment")
  
  # Rename columns based on the condition using rename_with()
  df <- df %>%
    rename_with(.cols = all_of(columns_to_rename),
                .fn = ~ paste0(condition_value, "_", .x))
  
  return(df)
})

# Function to perform inner join on chr, start, end, and strand
merge_dfs <- function(df1, df2) {
  inner_join(df1, df2, by = c("gene_name","feature_id","chr", "start", "end", "strand","pA_type"))
}

# Reduce the list of dataframes by performing inner joins
merged_df <- Reduce(merge_dfs, siCFIM25_files)

#identify control genes
control_genes = merged_df %>%
  dplyr::filter(abs(siCFIM25_DMSO_change_in_usage) < 0.1 & abs(siCFIM25_DMAi_change_in_usage) < 0.1 & siCFIM25_DMSO_mean_usage > 0.2 & siCFIM25_DMAi_mean_usage > 0.2)

#keep pairs only, i.e genes where both sites are unchanged
control_pairs = control_genes %>%
  group_by(gene_name) %>%
  dplyr::filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA'))

#filter to select only the genes that are significantly changed by siCFIM25 + DMSO
filtered_table = merged_df %>%
  dplyr::filter(
    (abs(siCFIM25_DMSO_change_in_usage) > 0.1 & siCFIM25_DMSO_twostep_feature_padj < 0.05))

#again only want to keep siCFIM25-shortened genes where both pA sites pass this threshold
filtered_pairs = filtered_table %>%
  group_by(gene_name) %>%
  dplyr::filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA'))  %>% # Ensure each gene group has at least one 'proximal pA' and one 'distal pA'
  dplyr::filter(siCFIM25_DMSO_prox_to_dist_shift < -0.2)

#categorise pA sites according to the ratio of proximal-to-distal shifts in siCFIM25+DMAi versus siCFIM25+DMSO
data_categorized <- filtered_pairs %>%
  mutate(
    category = case_when(
      # no mitigation
      (siCFIM25_DMAi_prox_to_dist_shift / siCFIM25_DMSO_prox_to_dist_shift) > 0.8 & 
      (siCFIM25_DMAi_prox_to_dist_shift / siCFIM25_DMSO_prox_to_dist_shift) < 1.2  ~ "no mitigation",
      
      #partial mitigation
      (siCFIM25_DMAi_prox_to_dist_shift / siCFIM25_DMSO_prox_to_dist_shift) < 0.6 & 
        (siCFIM25_DMAi_prox_to_dist_shift / siCFIM25_DMSO_prox_to_dist_shift) > 0.3  ~ "DMAi partial mitigation",
      
      #full mitigation
      (siCFIM25_DMAi_prox_to_dist_shift / siCFIM25_DMSO_prox_to_dist_shift) < 0.2
        ~ "DMAi full mitigation",
      
      #exacerbated
      (siCFIM25_DMAi_prox_to_dist_shift / siCFIM25_DMSO_prox_to_dist_shift) > 1.33  ~ "DMAi-exacerbated",
      
      # Default case if none of the above conditions are met
      TRUE ~ "undetermined"
    )
  )

#add UTR_type information so that we can focus solely on TUTR genes
directory <- "/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/TUTR_vs_ALE/sig_CSV/siCFIM25/"
siCFIM25_UTR_type_files <- list.files(directory, pattern = ".*all*", full.names = TRUE)
siCFIM25_UTR_type_list <- lapply(siCFIM25_UTR_type_files, read.csv)

#combine all the tables into one table
siCFIM25_sig_UTR_type_table = do.call(rbind, siCFIM25_UTR_type_list)

#do same with control pairs
directory <- "/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/TUTR_vs_ALE/ctrl_CSV/siCFIM25/"
siCFIM25_UTR_type_files <- list.files(directory, pattern = ".*all*", full.names = TRUE)
siCFIM25_UTR_type_table_list <- lapply(siCFIM25_UTR_type_files, read.csv)

#combine all the tables into one table
ctrl_UTR_type_table = do.call(rbind, siCFIM25_UTR_type_table_list)

#join consistent sites with UTR type tables and then only keep one row per pA site
data_categorized_with_UTR_type = data_categorized %>%
  left_join(siCFIM25_sig_UTR_type_table %>% dplyr::select(feature_id, UTR_type), by = "feature_id") %>%
  group_by(feature_id) %>% dplyr::slice(1) %>% ungroup()

ctrl_pairs_categorized_with_UTR_type = control_pairs %>%
  left_join(ctrl_UTR_type_table %>% dplyr::select(feature_id, UTR_type), by = "feature_id") %>%
  group_by(feature_id) %>% dplyr::slice(1) %>% ungroup() 

#filter to get TUTR genes
TUTR_data_categorized = data_categorized_with_UTR_type %>% dplyr::filter(UTR_type == 'TUTR') %>% group_by(gene_name) %>%
  dplyr::filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA'))

TUTR_control_pairs = ctrl_pairs_categorized_with_UTR_type %>% dplyr::filter(UTR_type == 'TUTR') %>% group_by(gene_name) %>%
  dplyr::filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA'))

# Make sure 'category' is a factor and data is sorted by category
TUTR_data_categorized_single <- TUTR_data_categorized %>%
  group_by(gene_name) %>% dplyr::slice(1) %>%
  dplyr::filter(category != 'undetermined') %>%
  mutate(category = factor(category, levels = c('no mitigation', 'DMAi partial mitigation', 'DMAi full mitigation', 'DMAi-exacerbated'))) %>%
  arrange(category)

# Prepare the annotation data frame
heatmap_data <- as.matrix(TUTR_data_categorized_single[, c("siCFIM25_DMSO_prox_to_dist_shift", "siCFIM25_DMAi_prox_to_dist_shift")])
annotations <- data.frame(Category = TUTR_data_categorized_single$category)
# Make sure the row names of 'heatmap_data' match the 'gene_name' for proper labeling
rownames(heatmap_data) <- TUTR_data_categorized_single$gene_name

# Prepare the annotations data frame
annotations <- data.frame(Category = TUTR_data_categorized_single$category)
rownames(annotations) <- TUTR_data_categorized_single$gene_name  # Ensure row names match those of heatmap_data
annotation_colors = list(
  Category = c("no mitigation" = "#469de1", "DMAi partial mitigation" = "#ff9e4e", "DMAi full mitigation" = "#ff594e", "DMAi-exacerbated" = '#075c75')  # Customize as needed
)

# Define the color palette and breaks, as before
max_val <- max(abs(heatmap_data), na.rm = TRUE)
colors <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-max_val, max_val, length.out = length(colors) + 1)

#define colour scheme
RdBu=colorRampPalette(rev(brewer.pal(9,'RdBu')))(30)  #extrapolate 9 colours palette out to 30
prox_to_dist_shift_breaks = seq(-0.7, 0.7, length.out = 31)  # Create breaks evenly spaced from -0.7 to 0.7

# Adjusting cell size and margins
cell_height <- 0.5  # Adjust to your need for better label visibility
cell_width <- 125   # Adjust based on your data size and display requirement
margins <- c(0, 25)  # Left and bottom margins

# Plot the heatmap with RdBu color scheme and row breaks
pheatmap(
  heatmap_data,
  clustering_distance_rows = 'binary', 
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = c("siCFIM25 + DMSO", "siCFIM25 + DMAi"),
  cellheight = cell_height,
  cellwidth = cell_width,
  color = RdBu,
  breaks = prox_to_dist_shift_breaks,
  fontsize_row = 1.25,
  annotation_row = annotations,
  annotation_colors = annotation_colors,
  angle_col = 45
)

#do siCFIM25 DMSO boxplot by category
data_for_siCFIM25_DMSO_boxplot = TUTR_data_categorized %>% filter(category != "undetermined" & category != "DMAi-exacerbated")

#factor the category variables
data_for_siCFIM25_DMSO_boxplot$category <- factor(
  data_for_siCFIM25_DMSO_boxplot$category,
  levels = c("no mitigation", "DMAi partial mitigation", "DMAi full mitigation")
)
  
ggplot(data_for_siCFIM25_DMSO_boxplot,
       aes(x = category,
           y = siCFIM25_DMSO_prox_to_dist_shift, fill = category)) +
  scale_fill_manual(values = c("no mitigation" = "#469de1","DMAi partial mitigation" = "#ff9e4e","DMAi full mitigation" = "#ff594e"))+
  geom_boxplot() +              # dots per point
  labs(
    x = "Category",
    y = "proximal to distal shift"
  ) +
  ylim(-1.25,0) +
  theme_classic(base_size = 14) +
  theme(
    text        = element_text(size = 18), 
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(
  siCFIM25_DMSO_prox_to_dist_shift ~ category,
  data = data_for_siCFIM25_DMSO_boxplot
)

dunn_test_result <- dunn.test(data_for_siCFIM25_DMSO_boxplot$siCFIM25_DMSO_prox_to_dist_shift, data_for_siCFIM25_DMSO_boxplot$category, method = "bonferroni") #full vs partial = 1.701669e-09, full vs no = 8.854442e-03, partial vs no = 1.530165e-03

#add difference in siCFIM25 DMAi - sicFIM25 DMSO proximal-to-distal shift (rather than ratio)
data_for_siCFIM25_DMSO_boxplot$DMAi_DMSO_difference = data_for_siCFIM25_DMSO_boxplot$siCFIM25_DMAi_prox_to_dist_shift - data_for_siCFIM25_DMSO_boxplot$siCFIM25_DMSO_prox_to_dist_shift
data_for_siCFIM25_DMSO_boxplot$DMSO_DMAi_difference = data_for_siCFIM25_DMSO_boxplot$siCFIM25_DMSO_prox_to_dist_shift - data_for_siCFIM25_DMSO_boxplot$siCFIM25_DMAi_prox_to_dist_shift

ggplot(data_for_siCFIM25_DMSO_boxplot,
       aes(x = category,
           y = DMAi_DMSO_difference, fill = category)) +
  scale_fill_manual(values = c("no mitigation" = "#469de1","DMAi partial mitigation" = "#ff9e4e","DMAi full mitigation" = "#ff594e"))+
  geom_boxplot() +              # dots per point
  labs(
    x = "Category",
    y = "siCFIM25+DMAi - siCFIM25+DMSO",
  ) +
  theme_classic(base_size = 14) +
  theme(
    text        = element_text(size = 18), 
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )

#isolate gene groups
no_mit_genes = TUTR_data_categorized %>% dplyr::filter(category == 'no mitigation')
DMAi_exac_genes = TUTR_data_categorized %>% dplyr::filter(category == 'DMAi-exacerbated')
DMAi_part_mit_genes = TUTR_data_categorized %>% dplyr::filter(category == 'DMAi partial mitigation')
DMAi_full_mit_genes = TUTR_data_categorized %>% dplyr::filter(category == 'DMAi full mitigation')
all_categories_genes = TUTR_data_categorized %>% dplyr::filter(category != 'undetermined')

#make list
gene_group_list = list(no_mit_genes,DMAi_exac_genes,DMAi_part_mit_genes,
                       DMAi_full_mit_genes,all_categories_genes,TUTR_control_pairs)

#add new columns defining window for fasta file
APA_sites_with_window_list = lapply(gene_group_list, function(df) {
    df <- df %>%
      mutate(
        region_start = if_else(strand == "+", end - 200, start - 200),
        region_end   = if_else(strand == "+", end + 200, start + 200)
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

#create dataframes to then convert into bed files
all_APA_sites_with_window_list = lapply(all_APA_sites_with_window_list,function(df){
  df = df[,c('chr','region_start','region_end','gene_name','feature_id','strand')]
  return(df)
})

#setwd and write to file
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/relaxed_thresholds/400")

#write bed
write.table(all_APA_sites_with_window_list[[1]], file = 'no_mit_distal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[2]], file = 'DMAi_exac_distal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[3]], file = 'DMAi_part_mit_distal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[4]], file = 'DMAi_full_mit_distal_sites.bed',
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[5]], file = 'siCFIM25_down_distal_sites.bed',
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(all_APA_sites_with_window_list[[7]], file = 'no_mit_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[8]], file = 'DMAi_exac_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[9]], file = 'DMAi_part_mit_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[10]], file = 'DMAi_full_mit_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[11]], file = 'siCFIM25_up_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##create expression matched groups
#read in deseq2 file
CFIM25_whole_gene_DE = read.csv("/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/DE_tables/siCFIM25_DMSO_whole_gene_deseq_results.csv")

#add basemean to gene group tables 
# Loop over each data frame in the list and perform left_join with CFIM25_whole_gene_DE
gene_group_with_basemean_list <- lapply(gene_group_list, function(df) {
  # Perform the left join to add Basemean column based on matching gene_name
  left_join(df, CFIM25_whole_gene_DE %>% dplyr::select(gene_name, baseMean), by = "gene_name")
})

# Calculate the percentile threshold of Basemean values in gene_group_list[[1]]
lower_threshold <- quantile(gene_group_with_basemean_list[[1]]$baseMean, 0.055, na.rm = TRUE)
upper_threshold_3 <- quantile(gene_group_with_basemean_list[[3]]$baseMean,0.91, na.rm = TRUE)
upper_threshold_4 <- quantile(gene_group_with_basemean_list[[4]]$baseMean, 0.965, na.rm = TRUE)

# Filter out rows where Basemean does not meet specified threshold
gene_group_with_basemean_list[[1]] <- gene_group_with_basemean_list[[1]] %>%
  dplyr::filter(baseMean >= lower_threshold | is.na(baseMean))

gene_group_with_basemean_list[[3]] <- gene_group_with_basemean_list[[3]] %>%
  dplyr::filter(baseMean <= upper_threshold_3 | is.na(baseMean))
gene_group_with_basemean_list[[4]] <- gene_group_with_basemean_list[[4]] %>%
  dplyr::filter(baseMean <= upper_threshold_4 | is.na(baseMean))

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

#create dataframes to then convert into bed files
all_APA_sites_with_window_list = lapply(all_APA_sites_with_window_list,function(df){
  df = df[,c('chr','region_start','region_end','gene_name','feature_id','strand')]
  return(df)
})

#setwd and write to file
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/expression_matched/2000/")

#write bed
write.table(all_APA_sites_with_window_list[[1]], file = 'no_mit_distal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[2]], file = 'DMAi_exac_distal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[3]], file = 'DMAi_part_mit_distal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[4]], file = 'DMAi_full_mit_distal_sites.bed',
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(all_APA_sites_with_window_list[[7]], file = 'no_mit_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[8]], file = 'DMAi_exac_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[9]], file = 'DMAi_part_mit_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_APA_sites_with_window_list[[10]], file = 'DMAi_full_mit_proximal_sites.bed', 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


### check proximal pA usage for expression matched groups
#only keep proximal pA rows and proximal pA usage + category columns
proximal_pA_usage_data_list = lapply(gene_group_with_basemean_list[1:4],function(df){
  df = df %>% dplyr::filter(pA_type == 'proximal pA') %>% dplyr::select(siCFIM25_DMSO_pA_usage_control,category) 
  return(df)
})

#combine into one dataframe for plotting
proximal_pA_usage_data = rbind(proximal_pA_usage_data_list[[1]],proximal_pA_usage_data_list[[3]],proximal_pA_usage_data_list[[4]])

#relevel category variable
proximal_pA_usage_data$category = factor(proximal_pA_usage_data$category, levels = c('no mitigation','DMAi partial mitigation','DMAi full mitigation'))

#ggplot
ggplot(proximal_pA_usage_data, aes(x = category, y = siCFIM25_DMSO_pA_usage_control, fill = category)) +
  geom_boxplot() +
  ylab("baseline proximal pA usage") +
  scale_fill_manual(values = c("no mitigation" = "#469de1","DMAi partial mitigation" = "#ff9e4e","DMAi full mitigation" = "#ff594e"))+
  ylim(0,0.9) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18),                        # Increase axis title size
    axis.text = element_text(size = 17),  
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.75), # Adjusted for readability
    legend.title = element_text(size = 16),                      # Increase legend title size
    legend.text = element_text(size = 14)                        # Increase legend text size
  ) 

#calculate mean value
no_mit_mean = mean(proximal_pA_usage_data$siCFIM25_DMSO_pA_usage_control[proximal_pA_usage_data$category == "no mitigation"],na.rm = TRUE) # 0.40
part_mit_mean = mean(proximal_pA_usage_data$siCFIM25_DMSO_pA_usage_control[proximal_pA_usage_data$category == "DMAi partial mitigation"], na.rm = TRUE) # 0.33
full_mit_mean = mean(proximal_pA_usage_data$siCFIM25_DMSO_pA_usage_control[proximal_pA_usage_data$category == "DMAi full mitigation"], na.rm = TRUE) # 0.34.

##get gene region bed files
#read in gencode file
gencode_gtf_path='/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/own_atlas/TUTR_vs_ALE/gencode.v45.annotation.gtf'

# Import the GTF file
gtf <- import(gencode_gtf_path)

# Filter for genes
genes <- gtf[gtf$type == "gene"]

# Extract metadata columns
genes_df <- as.data.frame(genes)

# Compute TSS based on strand
genes_df$TSS <- ifelse(genes_df$strand == "+", genes_df$start, genes_df$end)

# Define the function to merge a single data frame with genes_df
merge_gene_group <- function(df) {
  # Merge with genes_df on 'gene_name'
  merged_genes <- inner_join(df, genes_df, by = "gene_name")
  # Return the merged data frame
  return(merged_genes)
}

# Apply the merge function to each data frame in gene_group_list
merged_gene_groups <- lapply(gene_group_list, merge_gene_group)

add_distance_column <- function(df) {
  df %>%
    group_by(gene_name) %>%
    mutate(
      # Extract the start.x coordinate from the row where pA_type == 'proximal pA'
      prox_start = start.x[pA_type == 'proximal pA'][1],
      
      # Calculate TSS_to_pPA_distance using the prox_start value
      TSS_to_pPA_distance = if_else(
        strand.x == '+',
        prox_start - TSS,
        TSS - prox_start
      )
    ) %>%
    ungroup()
}

# Apply the function to each data frame in the list
merged_gene_groups_with_gene_distance <- lapply(merged_gene_groups, add_distance_column)

######make region bed files (for 10kb genes only)
#Only need to work with no mit, part mit, full mit
merged_gene_groups_with_gene_distance = merged_gene_groups_with_gene_distance[c(1,3,4)]

#filter to only keep genes > 10kb
merged_gene_groups_with_gene_distance_tenK = lapply(merged_gene_groups_with_gene_distance,function(df){
  df = df %>% filter(TSS_to_pPA_distance > 10000)
  return(df)
})

###first get 3' UTR region
merged_gene_groups_with_UTR_region <- lapply(
  merged_gene_groups_with_gene_distance_tenK,
  function(df) {
    df %>%
      group_by(gene_name) %>%
      mutate(
        UTR_region_start = if_else(
          strand.x == "+",
          min(end.x) - 1000,   # For plus-strand
          min(start.x) - 1000  # For minus-strand
        ),
        UTR_region_end = if_else(
          strand.x == "+",
          max(end.x) + 1000,   # For plus-strand
          max(start.x) + 1000  # For minus-strand
        )
      ) %>%
      ungroup()
  }
)

#only need one per gene
merged_gene_groups_with_UTR_region_single = lapply(merged_gene_groups_with_UTR_region,function(df){
  df = df %>%
    group_by(gene_name) %>%
    dplyr::slice(1) %>% ungroup()
})

# Define a function to prepare BED columns from a data frame
prepare_bed_df <- function(df) {
  bed_df <- data.frame(
    chr = df$seqnames,
    start = df$UTR_region_start , 
    end = df$UTR_region_end,
    name = df$gene_name,
    score = ".",
    strand = df$strand.x)  
  return(bed_df)
}

# Apply the function to each data frame in gene_groups_with_UTR_regions
UTR_bed_dfs <- lapply(merged_gene_groups_with_UTR_region_single, prepare_bed_df)

#### Define function to calculate upstream of gene regions

#only proximal row per gene
merged_gene_groups_with_gene_distance_tenK_proximal  = lapply(merged_gene_groups_with_gene_distance_tenK,function(df) {
  df = df %>% group_by(gene_name) %>% dplyr::filter(pA_type == 'proximal pA')
  return(df)
})

calculate_upstream_of_gene_regions <- function(df) {
  # Compute regions around TSS based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$TSS - 5000,   # For "+" strand, start is TSS - 5000 
    df$TSS + 2000             # For "-" strand, start is TSS + 2000
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$TSS - 2000,            # For "+" strand, end is TSS - 2000
    df$TSS + 5000             # For "-" strand, end is TSS + 5000
  )
  
  return(df)
}

# Define a function to prepare BED columns from a data frame
prepare_bed_df <- function(df) {
  bed_df <- data.frame(
    chr = df$seqnames,
    start = df$region_start , 
    end = df$region_end,
    name = df$gene_name,
    score = ".",
    strand = df$strand.x  
  )
  # Ensure 'start' is at least 0
  bed_df$start <- pmax(bed_df$start, 0)
  return(bed_df)
}

# Apply to each merged data frame
gene_groups_with_upstream_of_gene_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_proximal, calculate_upstream_of_gene_regions)

# Apply the function to each data frame in gene_groups_with_upstream_of_gene_regions
upstream_of_gene_bed_dfs <- lapply(gene_groups_with_upstream_of_gene_regions, prepare_bed_df)

#### Define function to calculate pre-TSS regions
calculate_pre_tss_regions <- function(df) {
  # Compute regions around TSS based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$TSS - 2000,   # For "+" strand, start is TSS - 2000 
    df$TSS + 1000             # For "-" strand, start is TSS + 1000
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$TSS - 1000,            # For "+" strand, end is TSS - 1000
    df$TSS + 2000             # For "-" strand, end is TSS + 2000
  )
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_pre_tss_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_proximal, calculate_pre_tss_regions)

# Apply the function to each data frame in gene_groups_with_tss_regions
pre_tss_bed_dfs <- lapply(gene_groups_with_pre_tss_regions, prepare_bed_df)

###### Define function to calculate TSS regions
calculate_tss_regions <- function(df) {
  # Calculate regions around TSS
  df$region_start <- pmax(df$TSS - 1000, 1)
  df$region_end <- df$TSS + 1000
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_tss_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_proximal, calculate_tss_regions)

# Apply the function to each data frame in gene_groups_with_tss_regions
tss_bed_dfs <- lapply(gene_groups_with_tss_regions, prepare_bed_df)

### Define function to calculate gene body regions
calculate_gene_body_regions <- function(df) {
  # Compute regions around TSS based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$TSS + 1000,  
    df$start.x + 1000          
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$end.x - 1000,      
    df$TSS - 1000         
  )
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_gene_body_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_proximal, calculate_gene_body_regions)

# Apply the function to each data frame in gene_groups_with_tss_regions
gene_body_bed_dfs <- lapply(gene_groups_with_gene_body_regions, prepare_bed_df)

##get distal pA sites for calculating downstream of gene regions
merged_gene_groups_with_gene_distance_tenK_distal  = lapply(merged_gene_groups_with_gene_distance_tenK,function(df) {
  df = df %>% group_by(gene_name) %>% dplyr::filter(pA_type == 'distal pA')
  return(df)
})

##get post distal pA site regions
calculate_post_dPA_regions <- function(df) {
  # Compute regions around pA based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$end.x + 1000,   
    df$start.x - 2000            
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$end.x + 2000,          
    df$start.x - 1000            
  )
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_post_dPA_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_distal, calculate_post_dPA_regions)

# Apply the function to each data frame in gene_groups_regions
post_dPA_bed_dfs <- lapply(gene_groups_with_post_dPA_regions, prepare_bed_df)

### Define function to calculate downstream of gene region
calculate_downstream_of_gene_regions <- function(df) {
  # Compute regions around pA based on strand
  df$region_start <- ifelse(
    df$strand.x == "+",
    df$end.x + 2000,   
    df$start.x - 5000            
  )
  
  df$region_end <- ifelse(
    df$strand.x == "+",
    df$end.x + 5000,          
    df$start.x - 2000            
  )
  
  return(df)
}

# Apply to each merged data frame
gene_groups_with_downstream_of_gene_regions <- lapply(merged_gene_groups_with_gene_distance_tenK_distal, calculate_downstream_of_gene_regions)

# Apply the function to each data frame in gene_groups_regions
downstream_of_gene_bed_dfs <- lapply(gene_groups_with_downstream_of_gene_regions, prepare_bed_df)

#combine all gene region bed files into a single list
gene_group_bed_files = c(upstream_of_gene_bed_dfs,pre_tss_bed_dfs,tss_bed_dfs,gene_body_bed_dfs,UTR_bed_dfs,post_dPA_bed_dfs,downstream_of_gene_bed_dfs)

#setwd to write bed files into 
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/boxplot_regions/")

#write bed
prefixes <- c('no_mit', 'part_mit', 'full_mit')
suffixes <- c('_upstream_of_gene.bed','_pre_TSS.bed', '_TSS.bed', '_gene_body.bed',
              '_UTR.bed', '_post_dPA.bed','_downstream_of_gene.bed')

# Loop over the data frames in bed_dfs
for (i in seq_along(gene_group_bed_files)) {
  # Get the current data frame
  bed_df <- gene_group_bed_files[[i]]
  
  # Determine the prefix index (1, 2, or 3)
  prefix_index <- ((i - 1) %% 3) + 1
  prefix <- prefixes[prefix_index]
  
  # Determine the suffix based on group number
  group_number <- ceiling(i / 3)
  suffix <- suffixes[group_number]
  
  # Construct the filename
  filename <- paste0(prefix, suffix)
  
  # Write the data frame to the file
  write.table(bed_df, file = filename, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}


