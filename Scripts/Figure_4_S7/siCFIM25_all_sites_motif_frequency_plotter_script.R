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

setwd("../../Data/Figure_3_S5_S6/fasta/1000nt_window/")

########################################################################
# Define Motifs
########################################################################

motifs <- c('TGTA')

########################################################################
# Function: Calculate Motif Running Average
########################################################################

calculate_motif_running_average <- function(sequences,
                                            motifs,
                                            short_name,
                                            window_size = 30) {
  
  n_seqs        <- length(sequences)
  L             <- Biostrings::width(sequences)[1]           # 1000 nt
  motif_patterns <- lapply(motifs, Biostrings::DNAString)
  motif_len <- nchar(motifs[[1]])   # 4 for UGUA-like 4-mers
  
  ## ------------------------------------------------------------------
  ## 1. presence/absence vector (how many sequences hit each position)
  ## ------------------------------------------------------------------
  presence_counts <- integer(L)                              # initialise to 0
  
  for (s in seq_len(n_seqs)) {
    hit_pos <- logical(L)                                    # FALSE everywhere
    for (pat in motif_patterns) {
      starts <- Biostrings::start(
        Biostrings::matchPattern(pat, sequences[[s]], fixed = TRUE))
      hit_pos[starts] <- TRUE                                # flag each start only once
    }
    presence_counts <- presence_counts + hit_pos             # add this sequence’s hits
  }
  
  ## keep only positions where a motif can start
  valid_L <- L - motif_len + 1
  presence_counts <- presence_counts[1:valid_L]
  
  ## ------------------------------------------------------------------
  ## 2. 30-nt running sum (left-aligned, step = 1 nt)
  ## ------------------------------------------------------------------
  kernel       <- rep(1L, window_size)
  running_sum  <- stats::filter(presence_counts, kernel, sides = 1)
  running_sum  <- running_sum[!is.na(running_sum)]           # drop leading NAs
  
  running_avg  <- running_sum / (window_size * n_seqs)
  
  positions <- seq_len(valid_L - window_size + 1) + floor(window_size/2) - 500      # –500 … +500
  
  data.frame(
    Position       = positions,
    motif_hits     = as.integer(running_sum),   # raw window hits
    RunningAverage = running_avg,
    site           = short_name,
    stringsAsFactors = FALSE
  )
}

########################################################################
# Function: Process Fasta Files
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
    # Tag each row with the 'site' name and pA_type
    result$site <- short_name
    result$pA_type <- ifelse(grepl("proximal", short_name), "proximal", "distal")
    
    # Store results and the count
    results_list[[file_name]] <- result
    counts_list[[short_name]] <- sequence_count
  }
  
  combined_results <- do.call(rbind, results_list)
  list(results = combined_results, counts = counts_list)
}

########################################################################
# Functions to List the Required Fasta Files
########################################################################

process_TUTR_directory_lengthened <- function(directory_name) {
  c(
    file.path(directory_name, paste0(directory_name, "_control_TUTR_distal_sites.fasta")),
    file.path(directory_name, paste0(directory_name, "_sig_up_TUTR_distal_sites.fasta")),
    file.path(directory_name, paste0(directory_name, "_control_TUTR_proximal_sites.fasta")),
    file.path(directory_name, paste0(directory_name, "_sig_down_TUTR_proximal_sites.fasta"))
  )
}

process_TUTR_directory_shortened <- function(directory_name) {
  c(
    file.path(directory_name, paste0(directory_name, "_control_TUTR_distal_sites.fasta")),
    file.path(directory_name, paste0(directory_name, "_sig_down_TUTR_distal_sites.fasta")),
    file.path(directory_name, paste0(directory_name, "_control_TUTR_proximal_sites.fasta")),
    file.path(directory_name, paste0(directory_name, "_sig_up_TUTR_proximal_sites.fasta"))
  )
}

########################################################################
# Example Usage: Gather Fasta File Paths
########################################################################

CFIM25_TUTR_fasta_files <- process_TUTR_directory_shortened("siCFIM25")

########################################################################
# Create Short Names for Plot Labels
########################################################################

CFIM25_short_names  <- c("control distal", "siCFIM25-shortened distal", 
                         "control proximal","siCFIM25-shortened proximal")

########################################################################
# Process the FASTA Files
########################################################################

CFIM25_TUTR_results = process_fasta_files(CFIM25_TUTR_fasta_files, CFIM25_short_names, motifs)

########################################################################
# Combine them into Lists
########################################################################

non_DMAi_TUTR_results_list <- list(CFIM25_TUTR_results)

########################################################################
# Reorder the pA_type Column (Proximal, Distal)
########################################################################

non_DMAi_TUTR_results_list <- lapply(non_DMAi_TUTR_results_list, function(data) {
  if ("results" %in% names(data)) {
    data$results$pA_type <- factor(data$results$pA_type, levels = c("proximal", "distal"),        
                                   labels = c("proximal pA", "distal pA"))
  }
  data
})

#CFIM25
siCFIM25_data <- non_DMAi_TUTR_results_list[[1]]  # i.e. CFIM25_TUTR_results
siCFIM25_df   <- siCFIM25_data$results
siCFIM25_counts <- siCFIM25_data$counts

# 1A) Create a 'condition' column that lumps proximal/distal for the legend
siCFIM25_df <- siCFIM25_df %>%
  mutate(
    condition = case_when(
      site %in% c("control distal", "control proximal") ~ "control",
      site %in% c("siCFIM25-shortened distal", "siCFIM25-shortened proximal") ~ "siCFIM25-shortened",
      TRUE ~ site
    )
  )

# 1B) pick only the "proximal" count
# so you don't double-sum. E.g., pick proximal if counts are the same.
picked_condition_counts <- c(
  "control"          = 0,
  "siCFIM25-shortened" = 0
)

for (nm in names(siCFIM25_counts)) {
  # We'll only pick the "proximal" entry
  # (Change to "distal" if you prefer.)
  if (grepl("proximal", nm)) {
    cond_key <- case_when(
      grepl("siCFIM25-shortened", nm) ~ "siCFIM25-shortened",
      grepl("control", nm)           ~ "control",
      TRUE                           ~ "other"
    )
    if (cond_key != "other") {
      picked_condition_counts[cond_key] <- siCFIM25_counts[[nm]]
    }
  }
}

# Build the final legend labels
labels_with_counts_cfim25 <- sapply(names(picked_condition_counts), function(cond) {
  paste0(cond, " (", picked_condition_counts[cond], ")")
})

siCFIM25_plot <- ggplot(siCFIM25_df, aes(x = Position, y = RunningAverage, color = condition)) +
  geom_line(aes(group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "siCFIM25",
    x = "Position Relative to Poly(A) Sites",
    y = "UGUA Motif Occurence",
    color = "Category"   # Legend title
  ) +
  lims(x = c(-180, 180), y = c(0.001, 0.018)) +
  theme_bw() +
  # Supply only 3 possible colors
  scale_color_manual(
    values = c(
      "control"           = "darkgrey",
      "siCFIM25-shortened" = "#469de1"
    ),
    labels = labels_with_counts_cfim25
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 18),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 14),
    legend.title  = element_text(size = 18),
    legend.text   = element_text(size = 13),
    strip.text    = element_text(size = 20)
  )

siCFIM25_plot
