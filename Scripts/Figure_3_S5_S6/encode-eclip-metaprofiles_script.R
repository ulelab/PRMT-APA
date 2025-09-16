#Author: Charlotte Capitanchik

#Load packages
library(data.table)
library(tidyverse)
library(ggthemes)
library(GenomicRanges)
library(genomation) #use version ‘1.19.1’, later version cause an error

# Function for reading in files and creating GRanges objects
read_and_prepare_bed <- function(file_path) {
  fread(file_path) %>%
    dplyr::rename(
      chr = V1, 
      start = V2,
      end = V3, 
      strand = V6
    ) %>%
    unique() %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

# List of DMAi vs non-DMAi transcript file paths
DMAi_vs_non_DMAi_bed_file_paths <- list(
  DMAi_proximal = "../../Data/Figure_3_S5_S6/bed/APA_3seq/expression_matched/2000nt_window/common_DMAi/common_DMAi_sig_down_TUTR_proximal_sites.bed",
  non_DMAi_proximal = "../../Data/Figure_3_S5_S6/bed/APA_3seq/expression_matched/2000nt_window/non_DMAi/non_DMAi_sig_down_TUTR_proximal_sites.bed",
  ctrl_proximal   = "../../Data/Figure_3_S5_S6/bed/APA_3seq/expression_matched/2000nt_window/common_DMAi/control_TUTR_proximal_sites.bed",
  DMAi_distal = "../../Data/Figure_3_S5_S6/bed/APA_3seq/expression_matched/2000nt_window/common_DMAi/common_DMAi_sig_up_TUTR_distal_sites.bed",
  non_DMAi_distal = "../../Data/Figure_3_S5_S6/bed/APA_3seq/expression_matched/2000nt_window/non_DMAi/non_DMAi_sig_up_TUTR_distal_sites.bed",
  ctrl_distal   = "../../Data/Figure_3_S5_S6/bed/APA_3seq/expression_matched/2000nt_window/common_DMAi/control_TUTR_distal_sites.bed"
)

# Read and prepare GRanges objects (change bed_file_path variable name to analyse different files)
gr_list <- lapply(DMAi_vs_non_DMAi_bed_file_paths, read_and_prepare_bed)

# Resize GRanges objects
resize_gr <- function(gr) {
  resize(gr, width = 400, fix = "center")
}

gr_resized_list <- lapply(gr_list, resize_gr)

# Load HepG2 or K562 eCLIP data file paths (switch file location and/or pattern to read in the different files)
eclip_beds <- list.files("/path/to/eclip/bed/xlink_files",
                         full.names = TRUE,
                         pattern = "*.bed.gz")
eclip_names <- basename(eclip_beds) %>% sub("-merged.xl.bed.gz","",.)

read_eclip_bed <- function(file) {
  readBed(file)
}

calculate_score_matrices <- function(eclip_gr, gr_resized_list, eclip_name) {
  # Initialize an empty list to store results
  score_matrices <- list()
  
  # Iterate over each condition
  for (condition_name in names(gr_resized_list)) {
    windows <- gr_resized_list[[condition_name]]
    
    # Calculate the score matrix
    sm <- ScoreMatrix(
      target = eclip_gr,
      windows = windows,
      strand.aware = TRUE,
      weight.col = NULL  # Adjust if you have a weight column
    )
    
    # Convert to data frame
    sm_df <- as.data.frame(sm)
    
    # Prepare data for plotting
    sm_long <- sm_df %>%
      pivot_longer(cols = everything(), names_to = "position", values_to = "value") %>%
      mutate(
        position = as.numeric(gsub("V", "", position)),
        condition = condition_name,
        sample = eclip_name
      ) %>%
      group_by(position, condition, sample) %>%
      summarize(mean = mean(value, na.rm = TRUE), .groups = "drop")
    
    # Store in the list
    score_matrices[[condition_name]] <- sm_long
  }
  
  # Combine all conditions
  combined_df <- bind_rows(score_matrices)
  return(combined_df)
}

# Initialize an empty list to store results
results_list <- list()

# Process each eCLIP file
for (i in seq_along(eclip_beds)) {
  eclip_file <- eclip_beds[i]
  eclip_name <- eclip_names[i]
  
  cat("Processing file:", eclip_name, "(", i, "of", length(eclip_beds), ")\n")
  
  # Read the eCLIP BED file
  eclip_gr <- read_eclip_bed(eclip_file)
  
  # Ensure the object is a GRanges
  if (!is(eclip_gr, "GRanges")) {
    cat("Skipping file:", eclip_name, "- not a GRanges object.\n")
    next
  }
  
  # Calculate score matrices
  sample_df <- calculate_score_matrices(eclip_gr, gr_resized_list, eclip_name)
  
  # Append to the results list
  results_list[[i]] <- sample_df
  
  # Clean up
  rm(eclip_gr, sample_df)
  gc()
}

# Combine all results into one data frame
final_to_plot <- bind_rows(results_list)

# Convert factors and set levels
final_to_plot$location <- ifelse(
  grepl("proximal", final_to_plot$condition), "proximal pA", "distal pA"
)

final_to_plot$location <- factor(final_to_plot$location, levels = c('proximal pA','distal pA'))

# Adjust condition factor levels
final_to_plot$condition <- gsub("_proximal|_distal", "", final_to_plot$condition)

# Adjust condition factor levels for DMAi vs non-DMAi data
final_to_plot$condition <- factor(final_to_plot$condition, levels = c("DMAi", "non_DMAi", "ctrl"))

#change condition to category
final_to_plot <- final_to_plot %>% dplyr::rename(category = condition)

## Count proximal sites per category (avoids double counting prox+distal)
cat_counts <- c(
  DMAi     = length(gr_list$DMAi_proximal),
  non_DMAi = length(gr_list$non_DMAi_proximal)
)

legend_breaks <- c("DMAi","non_DMAi")
legend_labels <- c(
  paste0("DMAi lengthened (",     cat_counts["DMAi"],     ")"),
  paste0("non-DMAi lengthened (", cat_counts["non_DMAi"], ")")
)

## Keep the factor order consistent
final_to_plot$category <- factor(final_to_plot$category, levels = legend_breaks)

#extract individual eCLIP sample data for plotting
final_to_plot_HepG2_CSTF2 = final_to_plot %>% dplyr::filter(sample == 'HepG2-CSTF2')
final_to_plot_HepG2_CSTF2T = final_to_plot %>% dplyr::filter(sample == 'HepG2-CSTF2T')
final_to_plot_K562_CSTF2T = final_to_plot %>% dplyr::filter(sample == 'K562-CSTF2T')
final_to_plot_K562_CPSF6 = final_to_plot %>% dplyr::filter(sample == 'K562-CPSF6')
final_to_plot_HepG2_FUBP3 = final_to_plot %>% dplyr::filter(sample == 'HepG2-FUBP3')

#remove ctrl individual eCLIP sample data for plotting
final_to_plot_HepG2_CSTF2_no_ctrl = final_to_plot_HepG2_CSTF2 %>% filter(category != 'ctrl')
final_to_plot_HepG2_CSTF2T_no_ctrl = final_to_plot_HepG2_CSTF2T %>% filter(category != 'ctrl')
final_to_plot_K562_CSTF2T_no_ctrl = final_to_plot_K562_CSTF2T %>% filter(category != 'ctrl')
final_to_plot_K562_CPSF6_no_ctrl = final_to_plot_K562_CPSF6 %>% filter(category != 'ctrl')
final_to_plot_HepG2_FUBP3_no_ctrl = final_to_plot_HepG2_FUBP3 %>% filter(category != 'ctrl')

ggplot(final_to_plot_K562_CPSF6_no_ctrl, aes(x = position, y = mean, color = category)) +
  stat_smooth(method = lm, formula = y ~ poly(x, 10), linewidth = 1.5) + 
  facet_wrap(~location, scales = "free_x", ncol = 2) +
  theme_bw() +
  geom_vline(xintercept = 200, linetype = "dotted") +
  scale_color_manual(
    values = c(DMAi = "#ff594e", non_DMAi = "#1A85FF"),
    breaks = legend_breaks,
    labels = legend_labels
  ) +
  scale_x_continuous(breaks = seq(0, 400, 100),
                     labels = seq(-200, 200, 100),
                     limits = c(0, 400)) +
  labs(x = "Position", y = "Mean", color = "Condition", title = "K562-CPSF6 eCLIP") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 20)
  )
