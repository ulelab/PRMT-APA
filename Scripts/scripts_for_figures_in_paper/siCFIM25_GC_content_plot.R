library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

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

## ── 1. locate FASTA files ──────────────────────────────────────────────
fasta_files <- list.files("fasta/1000", pattern = "\\.fasta$", full.names = TRUE)

## ── 2. GC-profile helper (40-nt sliding window) ────────────────────────
gc_profile <- function(seqs, window = 40) {
  L      <- width(seqs)[1]          
  nwin   <- L - window + 1           
  offset <- L/2                    
  
  nseq   <- length(seqs)
  gc_mat <- matrix(0L, nrow = nwin, ncol = nseq)
  
  for (j in seq_len(nseq)) {
    lf <- letterFrequencyInSlidingView(seqs[[j]],
                                       view.width = window,
                                       letters    = c("G", "C"))
    gc_mat[, j] <- rowSums(lf)
  }
  
  tibble(
    Position = seq_len(nwin) - offset,   # 1–961  →  –499…+461
    GC       = rowMeans(gc_mat) / window
  )
}

## ── 3. read each FASTA, compute profile, tack on meta-data ─────────────
all_df <- lapply(fasta_files, function(f) {
  fa     <- readDNAStringSet(f)
  prof   <- gc_profile(fa, window = 40)
  
  site   <- str_remove(basename(f), "\\.fasta$")
  tibble(
    site      = site,
    pA_type   = if_else(str_detect(site, "distal"), "distal pA", "proximal pA")
  ) %>% 
    dplyr::slice(rep(1, nrow(prof))) %>%         # repeat meta rows
    bind_cols(prof)
}) %>% bind_rows()

## ----------------------------------------------------------------------
##  3b.  assign each site to a Category and build the legend text
## ----------------------------------------------------------------------

# Check if the 'site' column exists in results_df
if ("site" %in% names(all_df)) {
  # Create the 'pA_type' column based on the content of 'site'
  all_df$pA_type <- ifelse(grepl("proximal", all_df$site), "proximal", "distal")
  
  # Convert 'pA_type' into a factor with a specified order
  all_df$pA_type = factor(
    all_df$pA_type,
    levels = c("proximal", "distal"),         # old levels
    labels = c("proximal pA", "distal pA")    # new labels
  )
}

#only keep no mit, part mit, and full mit sites
all_df = all_df %>%
  filter(site != 'DMAi_exac_distal_sites' & site != 'DMAi_exac_proximal_sites' & site != 'siCFIM25_down_distal_sites' & site != 'siCFIM25_up_proximal_sites')

plot_df <- all_df %>% filter(site != 'DMAi_exac_distal_sites') %>% filter(site != 'DMAi_exac_proximal_sites') %>%
  mutate(Category = case_when(
    site %in% c("DMAi_full_mit_distal_sites", "DMAi_full_mit_proximal_sites") ~ "DMAi full mitigation",
    site %in% c("DMAi_part_mit_distal_sites", "DMAi_part_mit_proximal_sites") ~ "DMAi part mitigation",
    site %in% c("no_mit_distal_sites", "no_mit_proximal_sites") ~ "no DMAi mitigation",
    TRUE ~ as.character(site)  # Default to existing site names if not in any group
  ))

##  counts per category   (proximal files only, as before)
cat_counts <- list(
  "DMAi full mitigation" = length(readDNAStringSet("fasta/1000/DMAi_full_mit_proximal_sites.fasta")),
  "DMAi part mitigation" = length(readDNAStringSet("fasta/1000/DMAi_part_mit_proximal_sites.fasta")),
  "no DMAi mitigation"   = length(readDNAStringSet("fasta/1000/no_mit_proximal_sites.fasta"))
)

labels_with_counts <- paste0(names(cat_counts), " (", unlist(cat_counts), ")")

plot_df$Category <- factor(plot_df$Category, levels = names(cat_counts))

## ----------------------------------------------------------------------
##  4.  GC plot with correct legend
## ----------------------------------------------------------------------
ggplot(plot_df,
       aes(x = Position, y = GC, colour = Category, group = site)) +
  geom_line(size = 1.8, alpha = 0.9) +
  facet_grid(~ pA_type, space = "free_y") +
  scale_colour_manual(
    values = c("DMAi full mitigation" = "#ff594e",
               "DMAi part mitigation" = "#ff9e4e",
               "no DMAi mitigation"   = "#469de1"),
    labels = labels_with_counts,
    breaks = names(cat_counts)
  ) +
  labs(x = "Position relative to poly(A) site",
       y = "Mean GC fraction (40-nt window)",
       colour = NULL) +
  lims(x = c(-270, 270), y = c(0.23, 0.5)) +
  theme_bw() +
  theme(
    axis.title   = element_text(size = 15),
    axis.text    = element_text(size = 13),
    legend.text  = element_text(size = 13),
    strip.text   = element_text(size = 18)
  )
