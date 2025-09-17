#### GC boxplot
#Define the directory containing the FASTA files
#setwd
setwd("../../Data/Figure_4_S7/fasta/")
fasta_dir <- "boxplot_regions/"

# Get the list of FASTA files
fasta_files <- list.files(path = fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

# Initialize an empty list to collect data frames
gc_content_list <- list()

# Define a mapping from file prefixes to condition labels
condition_map <- list(
  'no_mit' = 'no mit',
  'part_mit' = 'DMAi part mit',
  'full_mit' = 'DMAi full mit'
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

# Convert 'condition' and 'sequence_type' to factors for controlled plotting
gc_content_df$condition <- factor(gc_content_df$condition, levels = c('no mit', 'DMAi part mit', 'DMAi full mit'))
gc_content_df$sequence_type <- factor(gc_content_df$sequence_type, levels = c('upstream_of_gene','pre_TSS','TSS','gene_body','UTR','post_dPA','downstream_of_gene'))

# Create a new factor variable with the desired labels
gc_content_df <- gc_content_df %>%
  mutate(condition_label = factor(
    condition,
    levels = c('no mit', 'DMAi part mit', 'DMAi full mit'),
    labels = c('no mitigation', 'DMAi partial mitigation', 'DMAi full mitigation')
  ))

# Create boxplots with the new labels and specified colors
ggplot(gc_content_df, aes(x = sequence_type, y = gc_content, fill = condition_label)) +
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
      "no mitigation" = "#469de1",
      "DMAi partial mitigation" = "#ff9e4e",
      "DMAi full mitigation" = "#ff594e"
    )
  ) +
  scale_x_discrete(labels = c("upstream_of_gene" = "-5kb to -2kb","pre_TSS" = "pre-TSS (-2kb to -1kb)", "TSS" = "TSS",'gene_body' = 'Gene Body','UTR' = "3' UTR",'post_dPA' = "post 3'UTR (+1kb to +2kb)",'downstream_of_gene' = '+2kb to +5kb')) +
  guides(fill = guide_legend(title = NULL))

#do statistical test 
# keep just the two groups we want to compare
full_vs_none <- gc_content_df %>%
  filter(condition %in% c("DMAi full mit", "no mit"))

wilcox_tbl <- full_vs_none %>%                       # 9 rows → 9 tests
  group_by(sequence_type) %>%                        # one test per region
  summarise(
    p_raw  = wilcox.test(gc_content ~ condition)$p.value,
    n_full = sum(condition == "DMAi full mit"),
    n_none = sum(condition == "no mit"),
    med_full = median(gc_content[condition == "DMAi full mit"]),
    med_none = median(gc_content[condition == "no mit"]),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_raw, method = "BH"))      # FDR-correct

wilcox_tbl
