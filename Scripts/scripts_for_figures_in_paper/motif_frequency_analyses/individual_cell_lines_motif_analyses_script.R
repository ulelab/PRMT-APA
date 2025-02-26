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

setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/fasta/non_removal/1000/")

########################################################################
# Define Motifs
########################################################################

motifs <- c('TGTA')

########################################################################
# Function: Calculate Motif Running Average
########################################################################

calculate_motif_running_average <- function(sequences, motifs, short_name, window_size = 30) {
  num_sequences <- length(sequences)
  sequence_length <- nchar(sequences[[1]])
  positions <- seq(1, sequence_length - window_size + 1)
  presence <- integer(length(positions))
  
  # Convert motifs to DNAString objects
  motif_patterns <- lapply(motifs, DNAString)
  
  for (i in seq_along(positions)) {
    window_sequence <- subseq(sequences, start = i, end = i + window_size - 1)
    presence[i] <- sum(sapply(motif_patterns, function(motif) {
      sum(vcountPattern(motif, window_sequence) > 0)
    }))
  }
  
  running_average <- presence / window_size / num_sequences
  
  data.frame(
    Position = positions - 500,
    Presence = presence,
    RunningAverage = running_average,
    site = short_name
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

GP2D_TUTR_fasta_files <- process_TUTR_directory_lengthened("GP2D")
H838_TUTR_fasta_files = process_TUTR_directory_lengthened("H838")
HCTp53minus_TUTR_fasta_files = process_TUTR_directory_lengthened("HCT116_p53minus")
HCTp53plus_TUTR_fasta_files = process_TUTR_directory_lengthened("HCT116_p53plus")
LU99_TUTR_fasta_files = process_TUTR_directory_lengthened("LU99")
MCF7_TUTR_fasta_files = process_TUTR_directory_lengthened("MCF7")
PANC0403_TUTR_fasta_files = process_TUTR_directory_lengthened("PANC0403")
SKMEL2_TUTR_fasta_files = process_TUTR_directory_lengthened("SKMEL2")
SUM149PT_TUTR_fasta_files = process_TUTR_directory_lengthened("SUM149PT")
U2OS_TUTR_fasta_files = process_TUTR_directory_lengthened("U2OS")

CFIM25_TUTR_fasta_files <- process_TUTR_directory_shortened("siCFIM25")
PCF11_TUTR_fasta_files <- process_TUTR_directory_lengthened("siPCF11")
JTE607_TUTR_fasta_files <- process_TUTR_directory_lengthened("JTE607")

########################################################################
# Create Short Names for Plot Labels
########################################################################

DMAi_short_names    <- c("control distal", "DMAi-lengthened distal", 
                         "control proximal","DMAi-lengthened proximal")

CFIM25_short_names  <- c("control distal", "siCFIM25-shortened distal", 
                         "control proximal","siCFIM25-shortened proximal")

JTE607_short_names  <- c("control distal", "CPSF73i-lengthened distal",
                         "control proximal","CPSF73i-lengthened proximal")

PCF11_short_names   <- c("control distal", "siPCF11-lengthened distal", 
                         "control proximal","siPCF11-lengthened proximal")

########################################################################
# Process the FASTA Files
########################################################################

GP2D_TUTR_results <- process_fasta_files(GP2D_TUTR_fasta_files, DMAi_short_names, motifs)
H838_TUTR_results <- process_fasta_files(H838_TUTR_fasta_files, DMAi_short_names, motifs)
HCTp53minus_TUTR_results = process_fasta_files(HCTp53minus_TUTR_fasta_files, DMAi_short_names, motifs)
HCTp53plus_TUTR_results = process_fasta_files(HCTp53plus_TUTR_fasta_files, DMAi_short_names, motifs)
LU99_TUTR_results = process_fasta_files(LU99_TUTR_fasta_files, DMAi_short_names, motifs)
MCF7_TUTR_results = process_fasta_files(MCF7_TUTR_fasta_files, DMAi_short_names, motifs)
PANC0403_TUTR_results = process_fasta_files(PANC0403_TUTR_fasta_files, DMAi_short_names, motifs)
SUM149PT_TUTR_results = process_fasta_files(SUM149PT_TUTR_fasta_files, DMAi_short_names, motifs)
SKMEL2_TUTR_results = process_fasta_files(SKMEL2_TUTR_fasta_files, DMAi_short_names, motifs)
U2OS_TUTR_results = process_fasta_files(U2OS_TUTR_fasta_files, DMAi_short_names, motifs)

CFIM25_TUTR_results = process_fasta_files(CFIM25_TUTR_fasta_files, CFIM25_short_names, motifs)
PCF11_TUTR_results = process_fasta_files(PCF11_TUTR_fasta_files, PCF11_short_names, motifs)
JTE607_TUTR_results = process_fasta_files(JTE607_TUTR_fasta_files, JTE607_short_names, motifs)

########################################################################
# Combine them into Lists (DMAi vs non-DMAi)
########################################################################

DMAi_TUTR_results_list <-  list(GP2D_TUTR_results,H838_TUTR_results,HCTp53minus_TUTR_results,HCTp53plus_TUTR_results,
                                LU99_TUTR_results,MCF7_TUTR_results,PANC0403_TUTR_results,SUM149PT_TUTR_results,SKMEL2_TUTR_results,U2OS_TUTR_results)

non_DMAi_TUTR_results_list <- list(CFIM25_TUTR_results, PCF11_TUTR_results, JTE607_TUTR_results)

########################################################################
# Reorder the pA_type Column (Proximal, Distal)
########################################################################

DMAi_TUTR_results_list <- lapply(DMAi_TUTR_results_list, function(data) {
  if ("results" %in% names(data)) {
    data$results$pA_type <- factor(data$results$pA_type, levels = c("proximal", "distal"),        
                                   labels = c("proximal pA", "distal pA"))
  }
  data
})

non_DMAi_TUTR_results_list <- lapply(non_DMAi_TUTR_results_list, function(data) {
  if ("results" %in% names(data)) {
    data$results$pA_type <- factor(data$results$pA_type, levels = c("proximal", "distal"),        
                                   labels = c("proximal pA", "distal pA"))
  }
  data
})

# Example loop over DMAi_TUTR_results_list
processed_DMAi_list <- lapply(DMAi_TUTR_results_list, function(x) {
  # 'x' has two components:
  #   x$results = data frame of motif results
  #   x$counts  = named vector of counts
  
  # 1) Mutate the data frame to unify distal/proximal under a single 'condition' column
  df <- x$results %>%
    mutate(
      condition = case_when(
        site %in% c("control distal", "control proximal") ~ "control",
        site %in% c("DMAi-lengthened distal", "DMAi-lengthened proximal") ~ "DMAi-lengthened",
        TRUE ~ site  # fallback if something unexpected
      )
    )
  
  # 2) Build a named vector for the final counts, picking only proximal
  picked_condition_counts <- c(
    "control"         = 0,
    "DMAi-lengthened" = 0
  )
  
  cts <- x$counts  # This is your original counts list (named by site)
  
  # Loop through each named count
  for (nm in names(cts)) {
    # Only pick the "proximal" entry so you don't double-sum
    if (grepl("proximal", nm)) {
      cond_key <- case_when(
        grepl("DMAi-lengthened", nm) ~ "DMAi-lengthened",
        grepl("control", nm)         ~ "control",
        TRUE                         ~ "other"
      )
      if (cond_key != "other") {
        picked_condition_counts[cond_key] <- cts[[nm]]
      }
    }
  }
  
  # Build final legend labels
  labels_with_counts <- sapply(names(picked_condition_counts), function(cond) {
    paste0(cond, " (", picked_condition_counts[cond], ")")
  })
  
  # Return a small list with the updated data, plus new labels
  list(
    df                  = df,                   # updated results
    original_counts     = cts,                  # original counts if you need them
    labels_with_counts  = labels_with_counts    # final legend labels
  )
})

#plot
GP2D = ggplot(processed_DMAi_list[[1]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "GP2D",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[1]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )

H838 = ggplot(processed_DMAi_list[[2]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "H838",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[2]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )


HCTp53minus = ggplot(processed_DMAi_list[[3]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "HCT116 p53 -/-",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[3]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )


HCTp53plus = ggplot(processed_DMAi_list[[4]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "HCT116 p53 +/+",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[4]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )


LU99 = ggplot(processed_DMAi_list[[5]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "LU99",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[5]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )


MCF7 = ggplot(processed_DMAi_list[[6]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "MCF7",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[6]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )


PANC0403 = ggplot(processed_DMAi_list[[7]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "PANC0403",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[7]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )


SUM149PT = ggplot(processed_DMAi_list[[8]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "SUM149PT",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[8]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )

SKMEL2 = ggplot(processed_DMAi_list[[9]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "SK-MEL2",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[9]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )


U2OS = ggplot(processed_DMAi_list[[10]]$df, aes(x = Position, y = RunningAverage)) +
  geom_line(aes(color = condition, group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "U2OS",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"  # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.067)) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "control"         = "darkgrey",
      "DMAi-lengthened" = "#e53b3b"
    ),
    labels = processed_DMAi_list[[10]]$labels_with_counts
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 15)
  )

grid.arrange(GP2D,H838,HCTp53minus,HCTp53plus,LU99,MCF7,PANC0403,SKMEL2,SUM149PT,U2OS)


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
  lims(x = c(-180, 180), y = c(0.002, 0.0145)) +
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

#PCF11
siPCF11_data <- non_DMAi_TUTR_results_list[[2]]  # i.e. PCF11_TUTR_results
siPCF11_df   <- siPCF11_data$results
siPCF11_counts <- siPCF11_data$counts

# 1A) Create a 'condition' column that lumps proximal/distal for the legend
siPCF11_df <- siPCF11_df %>%
  mutate(
    condition = case_when(
      site %in% c("control distal", "control proximal") ~ "control",
      site %in% c("siPCF11-lengthened distal", "siPCF11-lengthened proximal") ~ "siPCF11-lengthened",
      TRUE ~ site
    )
  )

# 1B) pick only the "proximal" count
# so you don't double-sum. E.g., pick proximal if counts are the same.
picked_condition_counts <- c(
  "control"          = 0,
  "siPCF11-lengthened" = 0
)

for (nm in names(siPCF11_counts)) {
  # We'll only pick the "proximal" entry
  # (Change to "distal" if you prefer.)
  if (grepl("proximal", nm)) {
    cond_key <- case_when(
      grepl("siPCF11-lengthened", nm) ~ "siPCF11-lengthened",
      grepl("control", nm)           ~ "control",
      TRUE                           ~ "other"
    )
    if (cond_key != "other") {
      picked_condition_counts[cond_key] <- siPCF11_counts[[nm]]
    }
  }
}

# Build the final legend labels
labels_with_counts_siPCF11 <- sapply(names(picked_condition_counts), function(cond) {
  paste0(cond, " (", picked_condition_counts[cond], ")")
})

siPCF11_plot <- ggplot(siPCF11_df, aes(x = Position, y = RunningAverage, color = condition)) +
  geom_line(aes(group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "siPCF11",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"   # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.07)) +
  theme_bw() +
  # Supply only 3 possible colors
  scale_color_manual(
    values = c(
      "control"           = "darkgrey",
      "siPCF11-lengthened" = "#669b5f"
    ),
    labels = labels_with_counts_siPCF11
  ) +
  facet_grid(~pA_type, space = "free_y") +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 18),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 14),
    legend.title  = element_text(size = 18),
    legend.text   = element_text(size = 13),
    strip.text    = element_text(size = 18)
  )

siPCF11_plot

#JTE607
CPSF73i_data <- non_DMAi_TUTR_results_list[[3]]  # i.e. JTE607_TUTR_results
CPSF73i_df   <- CPSF73i_data$results
CPSF73i_counts <- CPSF73i_data$counts

# 1A) Create a 'condition' column that lumps proximal/distal for the legend
CPSF73i_df <- CPSF73i_df %>%
  mutate(
    condition = case_when(
      site %in% c("control distal", "control proximal") ~ "control",
      site %in% c("CPSF73i-lengthened distal", "CPSF73i-lengthened proximal") ~ "CPSF73i-lengthened",
      TRUE ~ site
    )
  )

# 1B) pick only the "proximal" count
# so you don't double-sum. E.g., pick proximal if counts are the same.
picked_condition_counts <- c(
  "control"          = 0,
  "CPSF73i-lengthened" = 0
)

for (nm in names(CPSF73i_counts)) {
  # We'll only pick the "proximal" entry
  # (Change to "distal" if you prefer.)
  if (grepl("proximal", nm)) {
    cond_key <- case_when(
      grepl("CPSF73i-lengthened", nm) ~ "CPSF73i-lengthened",
      grepl("control", nm)           ~ "control",
      TRUE                           ~ "other"
    )
    if (cond_key != "other") {
      picked_condition_counts[cond_key] <- CPSF73i_counts[[nm]]
    }
  }
}

# Build the final legend labels
labels_with_counts_CPSF73i <- sapply(names(picked_condition_counts), function(cond) {
  paste0(cond, " (", picked_condition_counts[cond], ")")
})

CPSF73i_plot <- ggplot(CPSF73i_df, aes(x = Position, y = RunningAverage, color = condition)) +
  geom_line(aes(group = site), linewidth = 1.3, alpha = 0.9) +
  labs(
    title = "CPSF73i",
    x = "Position Relative to Poly(A) Sites",
    y = "GC-rich Motif Occurence",
    color = "Category"   # Legend title
  ) +
  lims(x = c(-360, 360), y = c(0.01, 0.07)) +
  theme_bw() +
  # Supply only 3 possible colors
  scale_color_manual(
    values = c(
      "control"           = "darkgrey",
      "CPSF73i-lengthened" = "#669b5f"
    ),
    labels = labels_with_counts_CPSF73i
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

CPSF73i_plot

