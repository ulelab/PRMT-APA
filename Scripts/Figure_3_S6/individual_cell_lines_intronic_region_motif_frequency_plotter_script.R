library(Biostrings)
library(ggplot2)
library(tidyverse)
library(stringr)
library(ggh4x)

setwd("/Users/k2362866/Documents/AZ_postdoc/PRMT-APA_resubmission/Data/Figure_3_S6/fasta/1000nt_window")

# ============================================================
# Helper functions (unchanged from original)
# ============================================================

process_directory <- function(directory_name) {
  dir_path <- file.path(directory_name)
  fasta_files <- list.files(dir_path, pattern = "\\.fa$", full.names = TRUE)
  return(fasta_files)
}

calculate_motif_running_average <- function(sequences, motifs, short_name, window_size = 15) {
  n_seqs         <- length(sequences)
  L              <- Biostrings::width(sequences)[1]
  motif_patterns <- lapply(motifs, Biostrings::DNAString)
  motif_len      <- nchar(motifs[[1]])
  
  presence_counts <- integer(L)
  for (s in seq_len(n_seqs)) {
    hit_pos <- logical(L)
    for (pat in motif_patterns) {
      starts <- Biostrings::start(
        Biostrings::matchPattern(pat, sequences[[s]], fixed = TRUE))
      hit_pos[starts] <- TRUE
    }
    presence_counts <- presence_counts + hit_pos
  }
  
  valid_L         <- L - motif_len + 1
  presence_counts <- presence_counts[1:valid_L]
  
  running_sum <- stats::filter(presence_counts, rep(1L, window_size), sides = 2)
  running_sum <- running_sum[!is.na(running_sum)]
  running_avg <- running_sum / (window_size * n_seqs)
  positions   <- seq_len(valid_L - window_size + 1) + floor(window_size / 2) - 500
  
  data.frame(
    Position       = positions,
    motif_hits     = as.integer(running_sum),
    RunningAverage = running_avg,
    site           = short_name,
    stringsAsFactors = FALSE
  )
}

process_fasta_files <- function(fasta_files, motifs) {
  results_list <- list()
  counts_list  <- list()
  
  for (i in seq_along(fasta_files)) {
    fasta_sequences <- readDNAStringSet(fasta_files[i])
    file_name       <- basename(fasta_files[i])
    short_name      <- str_extract(file_name, "^.*(?=\\.fa$)")
    sequence_count  <- sum(grepl("^>", readLines(fasta_files[i]), perl = TRUE))
    result          <- calculate_motif_running_average(fasta_sequences, motifs, short_name)
    
    # ---- NEW: detect site_type from filename ----
    result$site_type <- case_when(
      grepl("intronic", file_name) ~ "intronic",
      grepl("UTR",      file_name) ~ "UTR",
      TRUE                         ~ "unknown"
    )
    
    results_list[[file_name]] <- result
    counts_list[[short_name]] <- sequence_count
  }
  
  combined_results <- do.call(rbind, results_list)
  return(list(results = combined_results, counts = counts_list))
}

# ============================================================
# Run
# ============================================================
motifs      <- c('GTGT','TGTG','GGTT','GTTG','TTGG','TGGT')
fasta_files <- process_directory("LU99_timecourse/")
fasta_output <- process_fasta_files(fasta_files, motifs)

results_cat <- fasta_output$results
counts_cat  <- fasta_output$counts

# ============================================================
# Annotate Category and site_type factor
# ============================================================
results_cat <- results_cat %>%
  mutate(
    Category = case_when(
      grepl("DMAi_increased",   site) ~ "DMAi increased",
      grepl("DMAi_unchanged",   site) ~ "DMAi unchanged",
      TRUE                            ~ site
    ),
    site_type = factor(
      site_type,
      levels = c("intronic", "UTR"),
      labels = c("intronic pA", "3'UTR pA")
    )
  )

# ============================================================
# Build counts lookup: one count per Category x site_type
# (use intronic counts as the n shown in legend, change to UTR if preferred)
# ============================================================
counts_df <- tibble(
  short_name = names(counts_cat),
  n          = unlist(counts_cat)
) %>%
  mutate(
    Category = case_when(
      grepl("DMAi_increased",   short_name) ~ "DMAi increased",
      grepl("DMAi_unchanged",   short_name) ~ "DMAi unchanged",
      TRUE                                  ~ short_name
    ),
    site_type = if_else(grepl("intronic", short_name), "intronic", "UTR")
  ) %>%
  filter(site_type == "intronic")   # use intronic n for legend labels

make_labels <- function(categories, counts_df) {
  map_chr(categories, function(cat) {
    n <- counts_df %>% filter(Category == cat) %>% pull(n)
    n <- if (length(n) == 0) "?" else n
    paste0(cat, " (n=", n, ")")
  })
}

# ============================================================
# Shared plot theme
# ============================================================
base_theme <- list(
  labs(x = "position relative to poly(A) site", y = "GU-rich Motif Occurrence"),
  lims(x = c(-40, 40), y = c(0.0025, 0.075)),
  theme_bw(),
  facet_grid(~ site_type, space = "free_y"),
  theme(
    axis.title  = element_text(size = 15),
    axis.text   = element_text(size = 14),
    legend.title = element_text(size = 17),
    legend.text  = element_text(size = 14),
    strip.text  = element_text(size = 18)
  )
)

# ============================================================
# Plot DMAi increased vs unchanged
# ============================================================
plot1_cats <- c(
  "DMAi increased",
  "DMAi unchanged"
)
plot1_colours <- c(
  "DMAi increased" = "#e05c4e",
  "DMAi unchanged"      = "#999999"
)

plot1_labels <- c(
  "DMAi increased" = paste0(
    "DMAi iAPA increased (",
    counts_df %>% filter(Category == "DMAi increased") %>% pull(n),
    ")"
  ),
  "DMAi unchanged" = paste0(
    "iAPA unchanged (",
    counts_df %>% filter(Category == "DMAi unchanged") %>% pull(n),
    ")"
  )
)
names(plot1_labels) <- plot1_cats

p1 <- results_cat %>%
  filter(Category %in% plot1_cats) %>%
  mutate(Category = factor(Category, levels = plot1_cats)) %>%
  ggplot(aes(x = Position, y = RunningAverage, color = Category)) +
  geom_line(linewidth = 2, alpha = 0.9) +
  scale_color_manual(
    values = plot1_colours,
    labels = plot1_labels
  ) +
  base_theme

print(p1)
