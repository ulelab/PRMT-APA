library(tidyverse)
library(data.table)
library(rstatix) 

#### GC boxplot
#Define the directory containing the FASTA files
#setwd
setwd("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines//Nobby_APA_analysis/common_atlas/resequenced/dedup/github/fasta/test_triple_removal/")
fasta_dir <- "boxplot_regions/DMAi_vs_non_DMAi/"

# Get the list of FASTA files
fasta_files <- list.files(path = fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

# Initialize an empty list to collect data frames
gc_content_list <- list()

# Define a mapping from file prefixes to condition labels
condition_map <- list(
  'control' = 'control',
  'DMAi' = 'DMAi lengthened',
  'non_DMAi' = 'non-DMAi lengthened'
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
      "DMAi lengthened" = "#ff594e",
      "non-DMAi lengthened" = "#1A85FF",
      "control" = "darkgrey"
    )
  ) +
  scale_x_discrete(labels = c("upstream_of_gene" = "-5kb to -2kb","pre_TSS" = "pre-TSS (-2kb to -1kb)", "TSS" = "TSS",'gene_body' = 'Gene Body','UTR' = "3' UTR",'post_dPA' = "post 3'UTR (+1kb to +2kb)",'downstream_of_gene' = '+2kb to +5kb')) +
  guides(fill = guide_legend(title = NULL))

#stats
gc_content_df <- gc_content_df %>%
  mutate(condition = factor(condition,
                            levels = c("control","DMAi lengthened","non-DMAi lengthened")))

# Per-region pairwise Wilcoxon comparing each group to control
pw_tbl <- gc_content_df %>%
  group_by(sequence_type) %>%
  pairwise_wilcox_test(
    gc_content ~ condition,
    ref.group = "control",      # compares control vs each other level
    p.adjust.method = "BH"      # FDR within each region
  ) %>%
  ungroup() %>%
  # keep just the two contrasts you care about
  filter(group2 %in% c("DMAi lengthened","non-DMAi lengthened")) %>%
  # (optional) FDR across *all* regions+comparisons:
  mutate(p.adj.global = p.adjust(p, method = "BH")) %>%
  mutate(comparison = paste(group1, "vs", group2))

     