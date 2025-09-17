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

# List of siCFIM25-transcript group (DMAi-sensitivity categorised) file paths, change path file to location of bed files you want to analyse
siCFIM25_bed_file_paths <- list(
  full_mit_proximal = "../../Data/Figure_4_S7/bed/APA_3seq/expression_matched/2000/DMAi_full_mit_proximal_sites.bed",
  part_mit_proximal = "../../Data/Figure_4_S7/bed/APA_3seq/expression_matched/2000/DMAi_part_mit_proximal_sites.bed",
  no_mit_proximal   = "../../Data/Figure_4_S7/bed/APA_3seq/expression_matched/2000/no_mit_proximal_sites.bed",
  full_mit_distal   = "../../Data/Figure_4_S7/bed/APA_3seq/expression_matched/2000/DMAi_full_mit_distal_sites.bed",
  part_mit_distal   = "../../Data/Figure_4_S7/bed/APA_3seq/expression_matched/2000/DMAi_part_mit_distal_sites.bed",
  no_mit_distal     = "../../Data/Figure_4_S7/bed/APA_3seq/expression_matched/2000/no_mit_distal_sites.bed"
)

# Read and prepare GRanges objects (change bed_file_path variable name to analyse different files)
gr_list <- lapply(siCFIM25_bed_file_paths, read_and_prepare_bed)

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

# Check lengths
length(eclip_beds)  
length(eclip_names)

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

# Adjust condition factor levels for siCFIM25 data
final_to_plot$condition <- factor(final_to_plot$condition, levels = c("no_mit", "part_mit", "full_mit"))

#change condition to category
final_to_plot <- final_to_plot %>% dplyr::rename(category = condition)

#extract individual eCLIP sample data for plotting
final_to_plot_HepG2_CSTF2 = final_to_plot %>% dplyr::filter(sample == 'HepG2-CSTF2')
final_to_plot_HepG2_CSTF2T = final_to_plot %>% dplyr::filter(sample == 'HepG2-CSTF2T')
final_to_plot_K562_CSTF2T = final_to_plot %>% dplyr::filter(sample == 'K562-CSTF2T')
final_to_plot_K562_CPSF6 = final_to_plot %>% dplyr::filter(sample == 'K562-CPSF6')
final_to_plot_HepG2_FUBP3 = final_to_plot %>% dplyr::filter(sample == 'HepG2-FUBP3')

# Plotting
ggplot(final_to_plot_HepG2_CSTF2, aes(x = position, y = mean, color = category)) +
  stat_smooth(
    aes(y = mean, x = position),
    method = lm,
    formula = y ~ poly(x, 10),
    linewidth = 1.2
  ) + 
  facet_wrap(~location, scales = "free_x",ncol = 2) +
  theme_bw() +
  geom_vline(xintercept = c(200), linetype = "dotted") +
  scale_color_manual(
    values = c(
      "no_mit" = "#469de1",
      "part_mit" = "#ff9e4e",
      "full_mit" = "#ff594e"
    ),
    labels = c(
      "DMAi No Mitigation",
      "DMAi Part Mitigation",
      "DMAi Full Mitigation"
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 400, by = 100),
    labels = seq(-200, 200, by = 100),
    limits = c(0, 400)
  ) +
  labs(
    x = "Position",
    y = "Mean",
    color = "Category",
    title = "HepG2-CSTF2 eCLIP"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 20)
  )

#filter by condition instead for plotting of all eCLIP samples
no_mit = final_to_plot %>% dplyr::filter(condition == 'no_mit')
part_mit = final_to_plot %>% dplyr::filter(condition == 'part_mit')
full_mit = final_to_plot %>% dplyr::filter(condition == 'full_mit')

#plot
ggplot(no_mit, aes(x = position, y = mean, color = sample)) +
  stat_smooth(
    aes(y = mean, x = position),
    method = lm,
    formula = y ~ poly(x, 10),
    linewidth = 1.5
  ) + 
  facet_wrap(~location, scales = "free_x",ncol = 2) +
  theme_bw() +
  geom_vline(xintercept = c(200), linetype = "dotted") +
  scale_x_continuous(
    breaks = seq(100, 300, by = 25),
    labels = seq(-100, 100, by = 25),
    limits = c(100, 300)
  ) +
  labs(
    x = "Position",
    y = "Mean",
    color = "Condition",
    title = "DMAi No Mitigation"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 20)
  )

