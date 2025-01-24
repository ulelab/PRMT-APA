library(Biostrings)
library(ggplot2)
library(tidyverse)
library(stringr)
library(ggh4x)

#setwd to location of fasta files directory
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/")

#function for reading in fasta files
process_directory <- function(directory_name) {
  # Define the path to the directory
  dir_path <- file.path(directory_name)  # Constructs the path to the directory
  
  # Get all files ending with .fasta in the specified directory
  fasta_files <- list.files(dir_path, pattern = "\\.fasta$", full.names = TRUE)
  
  return(fasta_files)
}

#execute fasta reading function
fasta_files = process_directory("fasta/1000/")

##### multiple motifs
# Define function to calculate running average of two motifs' presence/absence at each position
calculate_motif_running_average <- function(sequences, motifs, short_name, window_size = 10) {
  num_sequences <- length(sequences)
  sequence_length <- nchar(sequences[[1]])
  positions <- seq(1, sequence_length - window_size + 1)
  presence <- integer(length(positions))
  
  # Convert motifs to DNAString objects
  motif_patterns <- lapply(motifs, DNAString)
  
  for (i in seq_along(positions)) {
    # Sum the presence values for both motifs for the entire window
    window_sequence <- subseq(sequences, start = i, end = i + window_size - 1)
    presence[i] <- sum(sapply(motif_patterns, function(motif) {
      vcountPattern(motif, window_sequence) > 0
    }))
  }
  
  # Calculate the running average
  running_average = presence / window_size / num_sequences
  
  result <- data.frame(Position = positions - 500, Presence = presence, RunningAverage = running_average, site = short_name)
  return(result)
}

# Function to process fasta files and calculate motif presence
process_fasta_files <- function(fasta_files, motifs) {
  results_list <- list()
  
  for (i in seq_along(fasta_files)) {
    fasta_sequences <- readDNAStringSet(fasta_files[i])
    
    # Extract short_name using regex to capture the category and type
    file_name <- basename(fasta_files[i])
    short_name <- str_extract(file_name, "^.*(?=\\.fasta$)")  # Remove the .fasta extension
    
    result <- calculate_motif_running_average(fasta_sequences, motifs, short_name)
    
    # Add additional columns for plotting and analysis
    cat_type <- str_extract(file_name, ".*(?=_\\w+.fasta$)")  # Get the category like 'full_mit', etc.
    pA_type <- ifelse(grepl("proximal", file_name), "proximal", "distal")
    result$cat_type <- cat_type
    result$pA_type <- pA_type
    
    results_list[[short_name]] <- result
  }
  
  combined_results <- do.call(rbind, results_list)
  return(combined_results)
}

# Define motifs
motifs = c('AATAAA','ATTAAA')

#assess frequency of motifs 
results_cat = process_fasta_files(fasta_files, motifs)

# Check if the 'site' column exists in results_df
if ("site" %in% names(results_cat)) {
  # Create the 'pA_type' column based on the content of 'site'
  results_cat$pA_type <- ifelse(grepl("proximal", results_cat$site), "proximal", "distal")
  
  # Convert 'pA_type' into a factor with a specified order
  results_cat$pA_type <- factor(results_cat$pA_type, levels = c('proximal', 'distal'))  # Ensure 'proximal' comes first
}

#only keep no mit, part mit, and full mit sites
results_cat = results_cat %>%
  filter(site != 'DMAi_exac_distal_sites' & site != 'DMAi_exac_proximal_sites' & site != 'siCFIM25_down_distal_sites' & site != 'siCFIM25_up_proximal_sites')

#tidy up labels
results_cat <- results_cat %>%
  mutate(Category = case_when(
    site %in% c("DMAi_full_mit_distal_sites", "DMAi_full_mit_proximal_sites") ~ "DMAi full mitigation",
    site %in% c("DMAi_part_mit_distal_sites", "DMAi_part_mit_proximal_sites") ~ "DMAi part mitigation",
    site %in% c("no_mit_distal_sites", "no_mit_proximal_sites") ~ "no DMAi mitigation",
    TRUE ~ as.character(site)  # Default to existing site names if not in any group
  ))

##plot the graph
ggplot(results_cat, aes(x = Position, y = RunningAverage, color = Category)) +
  geom_line(linewidth = 2, alpha = 0.9) +
  labs(x = "Position", y = "Frequency") +
  lims(x = c(-50,20), y = c(0, 0.0725)) +
  theme_bw() +
  scale_color_manual(values = c("DMAi full mitigation" = "#ff594e",
                                "DMAi part mitigation" = "#ff9e4e",
                                "no DMAi mitigation" = "#469de1"),
                     labels = c("DMAi full mitigation" = "DMAi full mitigation",
                                "DMAi part mitigation" = "DMAi part mitigation",
                                "no DMAi mitigation" = "no DMAi mitigation")) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),  # Smaller title
    axis.title = element_text(size = 16),              # Smaller axis titles
    axis.text = element_text(size = 14),               # Smaller axis text
    legend.title = element_text(size = 18),            # Smaller legend titles
    legend.text = element_text(size = 14),  # Smaller legend text
    strip.text = element_text(size = 20)
  )

