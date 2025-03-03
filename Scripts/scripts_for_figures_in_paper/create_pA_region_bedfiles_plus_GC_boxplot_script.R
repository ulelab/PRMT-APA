library(tidyverse)
library(data.table)

##read in significant TUTR site files
#read in sig DMAi cell lines TUTR sites - change directory to the location of your TUTR significant site csv files are (generated from classification_of_APA_event_script.R script)
setwd("/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/panel_of_cancer_lines/APA_classified")
DMAi_data_frames = list.files(pattern = "^sig.*\\.csv")

# List to store data frames
sig_DMAi_TUTR_csv_files = list()

# Loop through each CSV file and read it
for (file in DMAi_data_frames) {
  df <- read_csv(file)
  sig_DMAi_TUTR_csv_files[[file]] <- df
}

#read in non-DMAi condition TUTR sites
setwd("/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig4/non_DMAi_datasets/APA_classified")
non_DMAi_data_frames = list.files(pattern = "^sig.*\\.csv")

# List to store data frames
sig_non_DMAi_TUTR_csv_files = list()

# Loop through each CSV file and read it
for (file in non_DMAi_data_frames) {
  df <- read_csv(file)
  sig_non_DMAi_TUTR_csv_files[[file]] <- df
}

#combine lists 
sig_TUTR_dfs = c(sig_DMAi_TUTR_csv_files,sig_non_DMAi_TUTR_csv_files)

## read in ctrl TUTR files
#read in ctrl DMAi cell lines TUTR sites - - change directory to the location of your TUTR ctrl site csv files are (generated from classification_of_APA_event_script.R script)
setwd("/ulelab/repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/panel_of_cancer_lines/APA_classified")
DMAi_data_frames = list.files(pattern = "^ctrl.*\\.csv")

# List to store data frames
ctrl_DMAi_TUTR_csv_files = list()

# Loop through each CSV file and read it
for (file in DMAi_data_frames) {
  df <- read_csv(file)
  ctrl_DMAi_TUTR_csv_files[[file]] <- df
}

#read in non-DMAi condition TUTR sites
setwd("/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig4/non_DMAi_datasets/APA_classified")
non_DMAi_data_frames = list.files(pattern = "^ctrl.*\\.csv")

# List to store data frames
ctrl_non_DMAi_TUTR_csv_files = list()

# Loop through each CSV file and read it
for (file in non_DMAi_data_frames) {
  df <- read_csv(file)
  ctrl_non_DMAi_TUTR_csv_files[[file]] <- df
}

#combine lists 
ctrl_TUTR_dfs = c(ctrl_DMAi_TUTR_csv_files,ctrl_non_DMAi_TUTR_csv_files)

#for motif analysis, only retain the most changing pA sites (and we want to keep only pairs that pass this threshold)
most_sig_TUTR_dfs = lapply(sig_TUTR_dfs,function(df){
  df = df %>% filter(abs(change_in_usage) > 0.15) %>%
    group_by(gene_name) %>%
    filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA')) %>%  # Ensure each gene group has at least one 'proximal pA' and one 'distal pA'
    mutate(sign_proximal = sign(change_in_usage[pA_type == 'proximal pA']),  # Determine the sign of avg_change_in_usage for proximal
           sign_distal = sign(change_in_usage[pA_type == 'distal pA'])) %>%  # Determine the sign of avg_change_in_usage for distal
    filter(sign_proximal != sign_distal) %>%  # Keep only those gene groups where the signs are different
    dplyr::select(-sign_proximal, -sign_distal) %>%  # Optionally remove the temporary sign columns
    ungroup() 
})

#for control sites, we want to make sure both sites in pair to be used 
both_used_ctrl_TUTR_dfs = lapply(ctrl_TUTR_dfs,function(df){
  df = df %>% filter(mean_usage > 0.15) %>% 
    group_by(gene_name) %>%
    filter(any(pA_type == 'proximal pA') & any(pA_type == 'distal pA'))  %>% ungroup()
  return(df)
})

## isolate distal and proximal
sig_TUTR_distal_sites_list = lapply(most_sig_TUTR_dfs,function(df){
  df = df %>% filter(pA_type == 'distal pA')
})

sig_TUTR_proximal_sites_list = lapply(most_sig_TUTR_dfs,function(df){
  df = df %>% filter(pA_type == 'proximal pA')
})

ctrl_TUTR_distal_sites_list = lapply(both_used_ctrl_TUTR_dfs,function(df){
  df = df %>% filter(pA_type == 'distal pA')
})

ctrl_TUTR_proximal_sites_list = lapply(both_used_ctrl_TUTR_dfs,function(df){
  df = df %>% filter(pA_type == 'proximal pA')
})

#find directionality of sig sites
up_TUTR_distal_sites_list = lapply(sig_TUTR_distal_sites_list,function(df){
  df = df %>% filter(change_in_usage > 0)
  return(df)
})

down_TUTR_distal_sites_list = lapply(sig_TUTR_distal_sites_list,function(df){
  df = df %>% filter(change_in_usage < 0)
  return(df)
})

up_TUTR_proximal_sites_list = lapply(sig_TUTR_proximal_sites_list,function(df){
  df = df %>% filter(change_in_usage > 0)
  return(df)
})

down_TUTR_proximal_sites_list = lapply(sig_TUTR_proximal_sites_list,function(df){
  df = df %>% filter(change_in_usage < 0)
  return(df)
})

#make list of common DMAi, non-DMAi, ctrl TUTR distal and proximal sites
all_TUTR_list = c(sig_TUTR_distal_sites_list,sig_TUTR_proximal_sites_list,up_TUTR_distal_sites_list,
                  down_TUTR_distal_sites_list,up_TUTR_proximal_sites_list,down_TUTR_proximal_sites_list,
                  ctrl_TUTR_distal_sites_list,ctrl_TUTR_proximal_sites_list)

### Define function to calculate pPA regions
calculate_PA_regions = function(df) {
  df <- df %>%
    mutate(
      region_start = if_else(strand == "+", end - 500, start - 500),
      region_end   = if_else(strand == "+", end + 500, start + 500)
    )
  return(df)
}

# Define a function to prepare BED columns from a data frame
prepare_bed_df <- function(df) {
  bed_df <- data.frame(
    chr = df$chr,
    start = df$region_start , 
    end = df$region_end,
    gene_name = df$gene_name,
    feature_id = df$feature_id,
    strand = df$strand  
  )
  # Ensure 'start' is at least 0
  bed_df$start <- pmax(bed_df$start, 0)
  return(bed_df)
}

# Apply to each merged data frame
all_TUTR_list_with_pA_regions <- lapply(all_TUTR_list, calculate_PA_regions)

# Apply the function to each data frame in gene_groups_with_tss_regions
all_TUTR_bed_files <- lapply(all_TUTR_list_with_pA_regions, prepare_bed_df)

# Function to write data frames to BED files with specified filenames and directories e.g.
setwd('/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/bed_files')

# Define the prefixes and suffixes
prefixes <- c('GP2D', 'H838', 'HCT116_p53minus', 'HCT116_p53plus', 'LU99', 'MCF7', 'PANC0403', 'SKMEL2', 'SUM149PT', 'U2OS', 'CBPi', 'CDK11i', 'CDK12i', 'CLKi_N', 'JTE607', 'METTL3i', 'siCFIM25', 'siELAVL1', 'siPCF11', 'siSam68', 'siSUPT4H1', 'siXRN2')

suffixes <- c("_sig_all_TUTR_distal_sites.bed", "_sig_all_TUTR_proximal_sites.bed", "_sig_up_TUTR_distal_sites.bed", "_sig_down_TUTR_distal_sites.bed", "_sig_up_TUTR_proximal_sites.bed", "_sig_down_TUTR_proximal_sites.bed", "_control_TUTR_distal_sites.bed", "_control_TUTR_proximal_sites.bed")

# Loop through each data frame in your list
for (i in 1:length(all_TUTR_bed_files)) {
  
  # Determine the suffix and prefix indices
  suffix_index <- ceiling(i / 22)
  prefix_index <- ((i - 1) %% 22) + 1
  
  # Get the corresponding prefix and suffix
  prefix <- prefixes[prefix_index]
  suffix <- suffixes[suffix_index]
  
  # Access the current data frame
  df <- all_TUTR_bed_files[[i]]
  
  # Create the directory if it doesn't exist
  dir.create(prefix, showWarnings = FALSE)
  
  # Construct the filename and filepath
  filename <- paste0(prefix, suffix)
  filepath <- file.path(prefix, filename)
  
  # Write the data frame to a bed file
  write.table(df, file = filepath, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

#### get gene region bed files for box plot
prox_gene_groups = c(down_TUTR_proximal_sites_list,ctrl_TUTR_proximal_sites_list)
distal_gene_groups = c(up_TUTR_distal_sites_list,ctrl_TUTR_distal_sites_list)
gene_groups <- Map(rbind, prox_gene_groups, distal_gene_groups)

#read in gencode file
gencode_gtf_path='gencode.v45.annotation.gtf'

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
merged_gene_groups <- lapply(gene_groups, merge_gene_group)

#add distance between TSS and prox pA
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
setwd("/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/bed_files/")

# Base prefixes (22)
base_prefixes <- c(
  "GP2D", "H838", "HCT116_p53minus", "HCT116_p53plus",
  "LU99", "MCF7", "PANC0403", "SKMEL2", "SUM149PT", "U2OS",
  "CBPi", "CDK11i", "CDK12i", "CLKi_N", "JTE607", "METTL3i",
  "siCFIM25", "siELAVL1", "siPCF11", "siSam68", "siSUPT4H1", "siXRN2"
)

# Suffixes (7)
suffixes <- c(
  "_upstream_of_gene.bed",
  "_pre_TSS.bed",
  "_TSS.bed",
  "_gene_body.bed",
  "_UTR.bed",
  "_post_dPA.bed",
  "_downstream_of_gene.bed"
)

# Total prefixes = 22 "base" + 22 "control" = 44
n_base <- length(base_prefixes)           # 22
n_total_prefixes <- n_base * 2           # 44
n_suffixes <- length(suffixes)           # 7

for (s in seq_len(n_suffixes)) {
  # Which suffix we are writing out
  current_suffix <- suffixes[s]
  
  # Loop over 44 prefix forms (the 22 base + 22 control)
  for (p in seq_len(n_total_prefixes)) {
    
    # Figure out the index in gene_group_bed_files
    df_index <- (s - 1) * n_total_prefixes + p
    bed_df   <- gene_group_bed_files[[df_index]]
    
    # Determine whether this is a "base" prefix or a "control" prefix
    base_index <- ((p - 1) %% n_base) + 1   # cycles 1..22
    is_control <- p > n_base               
    
    # The directory is always named after the base prefix
    dir_name   <- base_prefixes[base_index]
    
    # The file prefix is either "GP2D" or "GP2D_control", etc.
    if (is_control) {
      file_prefix <- paste0(dir_name, "_control")
    } else {
      file_prefix <- dir_name
    }
    
    # Create directory if it doesn't exist
    dir.create(dir_name, showWarnings = FALSE)
    
    # Build full path:  e.g. "GP2D/GP2D_control_downstream_of_gene.bed"
    filename <- file.path(dir_name, paste0(file_prefix, current_suffix))
    
    # Write out the bed file
    write.table(
      bed_df,
      file      = filename,
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
}

#### GC boxplot
#Define the directory containing the FASTA files - these need to be generated using getFasta from the bedtools suite, using the bed files created above us input
#setwd
setwd("fasta_files")
fasta_dir <- "boxplot_regions/JTE607" #e.g CPSF73i

# Get the list of FASTA files
fasta_files <- list.files(path = fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

# Initialize an empty list to collect data frames
gc_content_list <- list()

# Define a mapping from file prefixes to condition labels 
condition_map <- list(
  'JTE607_control' = 'control',
  'JTE607' = 'CPSF73i-lengthened'
)

# Loop over each FASTA file
for (fasta_file in fasta_files) {
  # Extract the base name of the file (without path)
  fasta_filename <- basename(fasta_file)
  
  # Remove the .fasta extension
  fasta_filename_no_ext <- sub("\\.fasta$", "", fasta_filename)
  
  # Extract the condition and sequence type using regular expressions
  if (grepl("_upstream_of_gene$", fasta_filename_no_ext)) {
    sequence_type <- 'upstream_of_gene'
    condition_key <- sub("_upstream_of_gene$", "", fasta_filename_no_ext)
  } else if (grepl("_downstream_of_gene$", fasta_filename_no_ext)) {
    sequence_type <- 'downstream_of_gene'
    condition_key <- sub("_downstream_of_gene$", "", fasta_filename_no_ext)
  } else if (grepl("_TSS_to_pPA$", fasta_filename_no_ext)) {
    sequence_type <- 'TSS_to_pPA'
    condition_key <- sub("_TSS_to_pPA$", "", fasta_filename_no_ext)
  } else if (grepl("_pPA_to_dPA$", fasta_filename_no_ext)) {
    sequence_type <- 'pPA_to_dPA'
    condition_key <- sub("_pPA_to_dPA$", "", fasta_filename_no_ext)
  } else if (grepl("_pre_TSS$", fasta_filename_no_ext)) {
    sequence_type <- 'pre_TSS'
    condition_key <- sub("_pre_TSS$", "", fasta_filename_no_ext)
  } else if (grepl("_TSS$", fasta_filename_no_ext)) {
    sequence_type <- 'TSS'
    condition_key <- sub("_TSS$", "", fasta_filename_no_ext)
  } else if (grepl("_gene_body$", fasta_filename_no_ext)) {
    sequence_type <- 'gene_body'
    condition_key <- sub("_gene_body$", "", fasta_filename_no_ext)
  } else if (grepl("_UTR$", fasta_filename_no_ext)) {
    sequence_type <- 'UTR'
    condition_key <- sub("_UTR$", "", fasta_filename_no_ext)
  } else if (grepl("_post_dPA$", fasta_filename_no_ext)) {
    sequence_type <- 'post_dPA'
    condition_key <- sub("_post_dPA$", "", fasta_filename_no_ext)
  } else {
    # Skip files that don't match expected patterns
    next
  }
  
  # Map the condition key to the condition label
  condition_label <- condition_map[[condition_key]]
  if (is.null(condition_label)) {
    condition_label <- condition_key  # Use original key if not in the map
  }
  
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file)
  
  # Calculate GC content for each sequence
  gc_counts <- letterFrequency(sequences, letters = c("G", "C"))
  sequence_lengths <- width(sequences)
  gc_content <- (rowSums(gc_counts) / sequence_lengths) * 100  # GC content as a proportion
  
  # Create a data frame
  df <- data.frame(
    sequence_name = names(sequences),
    gc_content = gc_content,
    condition = condition_label,
    sequence_type = sequence_type,
    stringsAsFactors = FALSE
  )
  
  # Append to the list
  gc_content_list[[length(gc_content_list) + 1]] <- df
}

# Combine all data frames in gc_content_list into one if needed
gc_content_df <- do.call(rbind, gc_content_list)

# Convert 'sequence_type' to factors for controlled plotting
gc_content_df$sequence_type <- factor(gc_content_df$sequence_type, levels = c('upstream_of_gene','pre_TSS','TSS','gene_body','UTR','post_dPA','downstream_of_gene'))

# Create boxplots with the new labels and specified colors
ggplot(gc_content_df, aes(x = sequence_type, y = gc_content, fill = condition)) +
  geom_boxplot() +
  ylab('GC Content (%)') +
  xlab('Gene region') +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 0.85,vjust = 0.9),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),  # Remove legend title
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  scale_fill_manual(
    values = c(
      "CPSF73i-lengthened" = "#669b5f",
      "control" = "darkgrey"
    )
  ) +
  scale_x_discrete(labels = c("upstream_of_gene" = "-5kb to -2kb","pre_TSS" = "pre-TSS (-2kb to -1kb)", "TSS" = "TSS",'gene_body' = 'Gene Body','UTR' = "3' UTR",'post_dPA' = "post 3'UTR (+1kb to +2kb)",'downstream_of_gene' = '+2kb to +5kb')) +
  guides(fill = guide_legend(title = NULL))

