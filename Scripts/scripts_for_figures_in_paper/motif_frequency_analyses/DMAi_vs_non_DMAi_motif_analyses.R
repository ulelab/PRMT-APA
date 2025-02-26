########################################################################
# Load Libraries
########################################################################

library(Biostrings)
library(ggplot2)
library(gridExtra)
library(ggridges)
library(tidyverse)
library(stringr)

########################################################################
# Set Working Directory
########################################################################

# Update this path to where your .fasta files are located
setwd("fasta_files/")

########################################################################
# Define Motifs
########################################################################

motifs <- c('AATAAA','ATTAAA')

########################################################################
# Function to Calculate Running Average of Motifs
########################################################################

calculate_motif_running_average <- function(sequences, motifs, short_name, window_size = 10) {
  num_sequences <- length(sequences)
  sequence_length <- nchar(sequences[[1]])
  positions <- seq(1, sequence_length - window_size + 1)
  presence <- integer(length(positions))
  
  # Convert motifs to DNAString objects
  motif_patterns <- lapply(motifs, DNAString)
  
  for (i in seq_along(positions)) {
    # Extract the window of sequences
    window_sequence <- subseq(sequences, start = i, end = i + window_size - 1)
    
    # Sum presence of all motifs within this window
    presence[i] <- sum(sapply(motif_patterns, function(motif) {
      sum(vcountPattern(motif, window_sequence) > 0)
    }))
  }
  
  # Calculate the running average
  running_average <- presence / window_size / num_sequences
  
  data.frame(
    Position = positions - 500,
    Presence = presence,
    RunningAverage = running_average,
    site = short_name
  )
}

########################################################################
# Function to Process Multiple Fasta Files
########################################################################

process_fasta_files <- function(fasta_files, short_names, motifs) {
  results_list <- list()
  counts_list <- list()
  
  for (i in seq_along(fasta_files)) {
    fasta_sequences <- readDNAStringSet(fasta_files[i])
    sequence_count <- sum(grepl("^>", readLines(fasta_files[i]), perl = TRUE))
    file_name <- basename(fasta_files[i])
    short_name <- short_names[i]
    
    result <- calculate_motif_running_average(fasta_sequences, motifs, short_name)
    
    # Tag the final data with the 'site' name and pA_type
    result$site <- short_name
    result$pA_type <- ifelse(grepl("proximal", short_name), "proximal", "distal")
    
    results_list[[file_name]] <- result
    counts_list[[short_name]] <- sequence_count
  }
  
  combined_results <- do.call(rbind, results_list)
  return(list(results = combined_results, counts = counts_list))
}

########################################################################
# Define File Retrieval Functions
########################################################################

process_common_TUTR_directory_lengthened <- function(directory_name) {
  c(
    paste0(directory_name, "/control_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_up_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/control_TUTR_proximal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_down_TUTR_proximal_sites.fasta")
  )
}

process_common_TUTR_directory_shortened <- function(directory_name) {
  c(
    paste0(directory_name, "/control_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_down_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/control_TUTR_proximal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_up_TUTR_proximal_sites.fasta")
  )
}

process_common_TUTR_directory_sig <- function(directory_name) {
  c(
    paste0(directory_name, "/control_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_all_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/control_TUTR_proximal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_all_TUTR_proximal_sites.fasta")
  )
}

process_non_DMAi_TUTR_directory_lengthened <- function(directory_name) {
  c(
    paste0(directory_name, "/", directory_name, "_sig_up_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_down_TUTR_proximal_sites.fasta")
  )
}

process_non_DMAi_TUTR_directory_shortened <- function(directory_name) {
  c(
    paste0(directory_name, "/", directory_name, "_sig_down_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_up_TUTR_proximal_sites.fasta")
  )
}

process_non_DMAi_TUTR_directory_sig <- function(directory_name) {
  c(
    paste0(directory_name, "/", directory_name, "_sig_all_TUTR_distal_sites.fasta"),
    paste0(directory_name, "/", directory_name, "_sig_all_TUTR_proximal_sites.fasta")
  )
}

########################################################################
# Define File Paths and Short Names
########################################################################

# Common DMAi short names
DMAi_short_names <- c("control distal", 
                      "DMAi-lengthened distal", 
                      "control proximal",
                      "DMAi-lengthened proximal")

# Non-DMAi short names
non_DMAi_short_names <- c("non-DMAi-lengthened distal", 
                          "non-DMAi-lengthened proximal")

# Gather the .fasta files
common_TUTR_fasta_files_lengthened <- process_common_TUTR_directory_lengthened("common_DMAi")
common_TUTR_fasta_files_shortened  <- process_common_TUTR_directory_lengthened("common_DMAi") 
common_TUTR_fasta_files_sig        <- process_common_TUTR_directory_sig("common_DMAi")

non_DMAi_TUTR_fasta_files_lengthened <- process_non_DMAi_TUTR_directory_lengthened("non_DMAi")
non_DMAi_TUTR_fasta_files_shortened  <- process_non_DMAi_TUTR_directory_lengthened("non_DMAi")
non_DMAi_TUTR_fasta_files_sig        <- process_non_DMAi_TUTR_directory_sig("non_DMAi")

########################################################################
# Process the Fasta Files
########################################################################

common_TUTR_results_lengthened  <- process_fasta_files(common_TUTR_fasta_files_lengthened,  DMAi_short_names, motifs)
common_TUTR_results_shortened   <- process_fasta_files(common_TUTR_fasta_files_shortened,   DMAi_short_names, motifs)
common_TUTR_results_sig         <- process_fasta_files(common_TUTR_fasta_files_sig,         DMAi_short_names, motifs)

non_DMAi_TUTR_results_lengthened <- process_fasta_files(non_DMAi_TUTR_fasta_files_lengthened, non_DMAi_short_names, motifs)
non_DMAi_TUTR_results_shortened  <- process_fasta_files(non_DMAi_TUTR_fasta_files_shortened,  non_DMAi_short_names, motifs)
non_DMAi_TUTR_results_sig        <- process_fasta_files(non_DMAi_TUTR_fasta_files_sig,        non_DMAi_short_names, motifs)

# Put results in a list for potential further use
TUTR_common_results_list <- list(
  common_TUTR_results_lengthened,
  common_TUTR_results_shortened,
  common_TUTR_results_sig,
  non_DMAi_TUTR_results_lengthened,
  non_DMAi_TUTR_results_shortened,
  non_DMAi_TUTR_results_sig
)

# Reorder proximal and distal levels if needed
TUTR_common_results_list <- lapply(TUTR_common_results_list, function(data) {
  data$results$pA_type <- factor(data$results$pA_type, levels = c("proximal", "distal"),        
                                 labels = c("proximal pA", "distal pA")    
  )
  data
})

########################################################################
# Combine Two Example Results for Plotting
########################################################################

# For demonstration, we'll combine just the "lengthened" sets
combined_results <- rbind(
  common_TUTR_results_lengthened$results,
  non_DMAi_TUTR_results_lengthened$results
)

combined_counts <- c(
  common_TUTR_results_lengthened$counts,
  non_DMAi_TUTR_results_lengthened$counts
)

# Create a single 'condition' column to unify proximal/distal
combined_results <- combined_results %>%
  mutate(
    # Group sites into 3 conditions:
    condition = case_when(
      site %in% c("DMAi-lengthened distal", "DMAi-lengthened proximal") ~ "DMAi-lengthened",
      site %in% c("non-DMAi-lengthened distal", "non-DMAi-lengthened proximal") ~ "non-DMAi-lengthened",
      site %in% c("control distal", "control proximal") ~ "control",
      TRUE ~ site
    )
  )

# Make sure pA_type is correct and ordered
combined_results$pA_type <- if_else(grepl("proximal", combined_results$site), "proximal", "distal")
combined_results$pA_type <- factor(combined_results$pA_type, levels = c("proximal", "distal"),        
                                   labels = c("proximal pA", "distal pA"))

########################################################################
# Compute Single Counts per Condition
########################################################################

condition_counts <- c(
  "DMAi-lengthened"       = 0,
  "non-DMAi-lengthened"   = 0,
  "control"                  = 0
)

for (nm in names(combined_counts)) {
  # Only use the 'proximal' count (change to 'distal' if you prefer)
  if (grepl("proximal", nm)) {
    cond_key <- case_when(
      grepl("non-DMAi-lengthened", nm) ~ "non-DMAi-lengthened",
      grepl("DMAi-lengthened", nm)     ~ "DMAi-lengthened",
      grepl("control", nm)                ~ "control",
      TRUE                             ~ "other"
    )
    condition_counts[cond_key] <- combined_counts[[nm]]
    }
  }


# Build the legend labels using these single (proximal) counts
labels_with_counts <- sapply(names(condition_counts), function(cond) {
  paste0(cond, " (", condition_counts[cond], ")")
})

########################################################################
# Plot
########################################################################

ggplot(combined_results, aes(x = Position, y = RunningAverage, color = condition)) +
  # group by 'site' so each site draws a separate line
  geom_line(aes(group = site), linewidth = 1.8, alpha = 0.8) +
  labs(
    x = "Position Relative to Poly(A) Sites",
    y = "AAUAAA/AUUAAA Motif Occurence",
    color = "Category"
  ) +
  lims(x = c(-40, 20), y = c(0, 0.075)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "DMAi-lengthened"       = "#e53b3b",
      "non-DMAi-lengthened"   = "darkgreen",
      "control"                  = "darkgrey"
    ),
    labels = labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14),
    axis.title       = element_text(size = 15),
    axis.text        = element_text(size = 14),
    legend.title     = element_text(size = 18),
    legend.text      = element_text(size = 14),
    strip.text       = element_text(size = 20)
  )
