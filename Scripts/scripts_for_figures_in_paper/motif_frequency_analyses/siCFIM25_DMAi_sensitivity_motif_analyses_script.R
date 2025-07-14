library(Biostrings)
library(ggplot2)
library(tidyverse)
library(stringr)
library(ggh4x)

#setwd to location of fasta files directory, e.g
setwd("path/to/fasta_files/")

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
  
  positions <- seq_len(valid_L - window_size + 1) - 500      # –500 … +500
  
  data.frame(
    Position       = positions,
    motif_hits     = as.integer(running_sum),   # raw window hits
    RunningAverage = running_avg,
    site           = short_name,
    stringsAsFactors = FALSE
  )
}


# Function to process fasta files and calculate motif presence
process_fasta_files <- function(fasta_files, motifs) {
  results_list <- list()
  counts_list <- list() 
  
  for (i in seq_along(fasta_files)) {
    fasta_sequences <- readDNAStringSet(fasta_files[i])
    file_name <- basename(fasta_files[i])
    
    # short_name could be the base name without .fasta
    short_name <- str_extract(file_name, "^.*(?=\\.fasta$)")
    
    # Count how many sequences are in this FASTA
    sequence_count <- sum(grepl("^>", readLines(fasta_files[i]), perl = TRUE))
    
    # Calculate motif running average (same as before)
    result <- calculate_motif_running_average(fasta_sequences, motifs, short_name)
    
    # Optionally store extra columns for plotting
    pA_type <- ifelse(grepl("proximal", file_name), "proximal", "distal")
    result$pA_type <- pA_type
    
    # Store the data frame in a list
    results_list[[file_name]] <- result
    
    # Store the counts in a named list
    # Key = short_name, Value = sequence_count
    counts_list[[short_name]] <- sequence_count
  }
  
  # Combine results into one data frame
  combined_results <- do.call(rbind, results_list)
  
  # Return BOTH the combined data frame and the named counts list
  return(list(results = combined_results, counts = counts_list))
}

# Define motifs
motifs = c('GGTT','TTGG','GTGT','TGTG','GTTG','TGGT')

#assess frequency of motifs 
fasta_output <- process_fasta_files(fasta_files, motifs)

results_cat <- fasta_output$results  # big data frame of all sites
counts_cat  <- fasta_output$counts   # named list of sequence counts

# Check if the 'site' column exists in results_df
if ("site" %in% names(results_cat)) {
  # Create the 'pA_type' column based on the content of 'site'
  results_cat$pA_type <- ifelse(grepl("proximal", results_cat$site), "proximal", "distal")
  
  # Convert 'pA_type' into a factor with a specified order
  results_cat$pA_type = factor(
    results_cat$pA_type,
    levels = c("proximal", "distal"),         # old levels
    labels = c("proximal pA", "distal pA")    # new labels
  )
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

# Compute Single Counts per Condition
condition_counts <- c(
  "DMAi full mitigation"       = 0,
  "DMAi part mitigation"   = 0,
  "no DMAi mitigation"                  = 0
)

for (nm in names(counts_cat)) {
  # Only use the 'proximal' count (change to 'distal' if you prefer)
  if (grepl("proximal", nm)) {
    cond_key <- case_when(
      grepl("DMAi_full_mit", nm) ~ "DMAi full mitigation",
      grepl("DMAi_part_mit", nm)     ~ "DMAi part mitigation",
      grepl("no_mit", nm)                ~ "no DMAi mitigation",
      TRUE                             ~ "other"
    )
    condition_counts[cond_key] <- counts_cat[[nm]]
  }
}

# Build the legend labels using these single (proximal) counts
labels_with_counts <- sapply(names(condition_counts), function(cond) {
  paste0(cond, " (", condition_counts[cond], ")")
})

##plot the graph
ggplot(results_cat, aes(x = Position, y = RunningAverage, color = Category)) +
  geom_line(linewidth = 2, alpha = 0.9) +
  labs(x = "Position Relative to Poly(A) Sites", y = "GU-rich Motif Occurence") +
  lims(x = c(-40, 40), y = c(0.005, 0.092)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "DMAi full mitigation" = "#ff594e",
      "DMAi part mitigation" = "#ff9e4e",
      "no DMAi mitigation"   = "#469de1"
    ),
    labels = labels_with_counts
  ) +
  # Facet with new factor labels
  facet_grid(~ pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 13),
    axis.title   = element_text(size = 15),
    axis.text    = element_text(size = 14),
    legend.title = element_text(size = 17),
    legend.text  = element_text(size = 14),
    strip.text   = element_text(size = 18)
  )
