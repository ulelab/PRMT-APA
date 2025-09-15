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
setwd("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/fasta/test_triple_removal/1000/")

########################################################################
# Define Motifs
########################################################################

motifs <- c('GGTT','TTGG','GTTG','TGGT','GTGT','TGTG')

########################################################################
# Function to Calculate Running Average of Motifs
########################################################################

##############################################################################
#  calculate_motif_running_average 
##############################################################################
calculate_motif_running_average <- function(sequences,
                                            motifs,
                                            short_name,
                                            window_size = 15) {
  
  n_seqs        <- length(sequences)
  L             <- Biostrings::width(sequences)[1]           
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
  geom_line(aes(group = site), linewidth = 1.5, alpha = 0.8) +
  labs(
    x = "Position Relative to Poly(A) Sites",
    y = "GU-rich motif occurence",
    color = "Category"
  ) +
  lims(x = c(-40, 40), y = c(0.005, 0.0925)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "DMAi-lengthened"       = "#ff594e",
      "non-DMAi-lengthened"   = "#1A85FF",
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
    strip.text       = element_text(size = 18)
  )

############################################################
##  Significance test 
############################################################

## ---- 1.  define window & build 2×2 table  ------------------------------
## ── 1.  Pull the counts for the window that starts at position 4 ─────────
start_pos   <- 4        # first nt of the 15-nt window

# Presence counts (k) in each condition
test_df <- combined_results %>%
  filter(Position == start_pos,                 
         pA_type == "distal pA",               # distal only
         condition %in% c("DMAi-lengthened",
                          "non-DMAi-lengthened"))

# total motif occurrences in that window, per group
k_full <- test_df$motif_hits[test_df$condition == "DMAi-lengthened"]
k_none <- test_df$motif_hits[test_df$condition == "non-DMAi-lengthened"]

# number of sequences (i.e. exposures) in each distal FASTA
n_full <- combined_counts[["DMAi-lengthened distal"]]
n_none <- combined_counts[["non-DMAi-lengthened distal"]]

## ── 2. Poisson rate test ───────────────────────────────────────────────
poisson_res <- poisson.test(
  x       = c(k_full, k_none),    # observed motif counts
  T       = c(n_full, n_none),    # exposure = # of sequences
  alternative = "two.sided"
)
print(poisson_res) #p-value = 6.191e-06

