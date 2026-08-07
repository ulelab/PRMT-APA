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
setwd("Data/Figure_4_S7/fasta/1000nt_window/")

########################################################################
# Define Motifs
########################################################################

motifs <- c('GGTT','GTGT','TTGG','TGTG','GTTG','TGGT')

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

fasta_files <- c(
  "PCF11_not_mitigated_distal_sites.fa",
  "PCF11_not_mitigated_proximal_sites.fa",
  "PCF11_strong_mitigation_distal_sites.fa",
  "PCF11_strong_mitigation_proximal_sites.fa"
)

short_names <- c(
  "not mitigated distal",
  "not mitigated proximal",
  "strong mitigation distal",
  "strong mitigation proximal"
)

results_obj <- process_fasta_files(
  fasta_files,
  short_names,
  motifs
)


combined_results <- results_obj$results
combined_counts  <- results_obj$counts

combined_results <- combined_results %>%
  mutate(
    condition = case_when(
      site %in% c(
        "not mitigated distal",
        "not mitigated proximal"
      ) ~ "not mitigated",
      
      site %in% c(
        "strong mitigation distal",
        "strong mitigation proximal"
      ) ~ "strong mitigation"
    )
  )

combined_results$pA_type <- if_else(
  grepl("proximal", combined_results$site),
  "proximal",
  "distal"
)

combined_results$pA_type <- factor(
  combined_results$pA_type,
  levels = c("proximal","distal"),
  labels = c("proximal pA","distal pA")
)

condition_counts <- c(
  "not mitigated" = 0,
  "strong mitigation" = 0
)

for (nm in names(combined_counts)) {
  
  if (grepl("proximal", nm)) {
    
    cond_key <- case_when(
      grepl("not mitigated", nm) ~ "not mitigated",
      grepl("strong mitigation", nm) ~ "strong mitigation"
    )
    
    condition_counts[cond_key] <- combined_counts[[nm]]
  }
}

labels_with_counts <- sapply(
  names(condition_counts),
  function(cond) {
    paste0(
      cond,
      " (",
      condition_counts[cond],
      ")"
    )
  }
)

ggplot(
  combined_results,
  aes(
    x = Position,
    y = RunningAverage,
    color = condition
  )
) +
  geom_line(
    aes(group = site),
    linewidth = 1.5,
    alpha = 0.8
  ) +
  labs(
    x = "Position Relative to Poly(A) Site",
    y = "GU-rich motif occurrence",
    color = "Category"
  ) +
  lims(
    x = c(-40,40)
  ) +
  scale_color_manual(
    values = c(
      "not mitigated" = "#d73027",
      "strong mitigation" = "#4575b4"
    ),
    labels = c(
      paste0("Not PCF11-oe mitigated (", condition_counts["not mitigated"], ")"),
      paste0("Strong PCF11-oe mitigation (", condition_counts["strong mitigation"], ")")
    )
  ) +
  facet_grid(~pA_type) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

