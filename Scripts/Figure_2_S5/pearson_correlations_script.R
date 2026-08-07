library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)

# Set the working directory to the directory containing the CSV files
setwd("Data/Figure_2_S5/APA_CSV_files/panel_of_cancer_lines/positional_info")

#read in files
csv_files <- list.files(pattern = "top2*")

# Read each file, keep a STANDARDISED column name, key the list by cell line
DMAi_files <- lapply(csv_files, function(file) {
  df <- read.csv(file) %>%
    dplyr::select(gene_name, feature_id, change_in_usage, condition, prox_to_dist_shift)
  df
})

# name each element by its condition (cell line), not by filename
line_names <- vapply(DMAi_files, function(df) unique(df$condition)[1], character(1))
names(DMAi_files) <- sub("^top2_", "", line_names)

line_names[line_names == "DMAi_72hrs"] <- "RPE-1"
names(DMAi_files) <- line_names

# and fix the condition column inside that df too, for consistency
DMAi_files[["RPE-1"]]$condition <- "RPE-1"

# join two lines by feature, return a tidy 2-column table
get_pair <- function(line_x, line_y, files = DMAi_files) {
  x <- files[[line_x]] %>% dplyr::select(gene_name, feature_id, x = change_in_usage)
  y <- files[[line_y]] %>% dplyr::select(feature_id, y = change_in_usage)
  inner_join(x, y, by = "feature_id")
}

plot_pair <- function(line_x, line_y, files = DMAi_files, lims = c(-0.6, 0.6)) {
  d <- get_pair(line_x, line_y, files)
  r <- cor(d$x, d$y, use = "complete.obs", method = "pearson")
  
  ggplot(d, aes(x = x, y = y)) +
    geom_hex(bins = 90, aes(fill = after_stat(count))) +
    scale_fill_viridis_c(option = "C") +
    geom_smooth(method = "lm", colour = "red", se = TRUE) +
    labs(x = paste0("Change in PAS usage: ", line_x),
         y = paste0("Change in PAS usage: ", line_y),
         title = paste0("DMAi APA - ", line_y, " vs ", line_x)) +
    coord_equal(xlim = lims, ylim = lims) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 26),
          axis.title = element_text(size = 18),
          axis.text  = element_text(size = 22),
          legend.position = "right") +
    annotate("text", x = lims[1] * 0.7, y = lims[2] * 0.75,
             label = sprintf("r = %.2f", r),
             fontface = "bold", size = 13, colour = "red")
}

# now every comparison is one line, by NAME:
plot_pair("LU99", "MCF7")
plot_pair("PANC0403", "GP2D")
plot_pair("H838", "SUM149PT")
plot_pair("HCT116_p53plus", "U2OS")

# names of the cancer lines (everything except RPE-1)
cancer_lines <- setdiff(names(DMAi_files), "RPE-1")

# stack all cancer lines into one long table
cancer_long <- purrr::map_dfr(cancer_lines, function(ln) {
  DMAi_files[[ln]] %>%
    dplyr::select(gene_name, feature_id, change_in_usage) %>%
    dplyr::mutate(line = ln)
})

# per-feature mean across whatever lines have that feature (na.rm handles absence)
cancer_mean <- cancer_long %>%
  dplyr::group_by(gene_name, feature_id) %>%
  dplyr::summarise(
    mean_change = mean(change_in_usage, na.rm = TRUE),
    n_lines = dplyr::n(),           # how many lines contributed
    .groups = "drop"
  )

rpe1_vs_mean <- cancer_mean %>%
  inner_join(
    DMAi_files[["RPE-1"]] %>% dplyr::select(feature_id, rpe1_change = change_in_usage),
    by = "feature_id"
  )

r_mean <- cor(rpe1_vs_mean$mean_change, rpe1_vs_mean$rpe1_change, use = "complete.obs")

ggplot(rpe1_vs_mean, aes(x = mean_change, y = rpe1_change)) +
  geom_hex(bins = 90, aes(fill = after_stat(count))) +
  scale_fill_viridis_c(option = "C") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40", linewidth = 0.8) +
  geom_smooth(method = "lm", colour = "red", se = TRUE) +
  labs(x = "Mean change in PAS usage: cancer lines",
       y = "Change in PAS usage: RPE-1",
       title = "DMAi APA - RPE-1 vs mean of cancer lines") +
  coord_equal(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 26),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text  = element_text(size = 15)) +
  annotate("text", x = -0.35, y = 0.45, label = sprintf("r = %.2f", r_mean),
           fontface = "bold", size = 10, colour = "red")
           
### lung cell line vs organoid comparison
#setwd for organoids
setwd("Data/Figure_2_S5/APA_CSV_files/patient_derived_organoids/positional_info/")
csv_files <- list.files(pattern = "\\.csv$")
csv_files <- csv_files[grepl("top2", csv_files) & grepl("DMAi", csv_files)]

# Read each CSV file and store them in a list called 'DMAi_files'
organoid_files <- lapply(csv_files, function(file) {
  df <- read.csv(file)  # Or use read_csv(file) for faster loading with readr
  
  # Select necessary columns
  df <- df %>%
    dplyr::select(gene_name, feature_id, change_in_usage, condition, prox_to_dist_shift)
  
  # Assume 'condition' is the same across all rows, take the first one
  condition_value <- unique(df$condition)[1]
  
  # Rename 'change_in_usage' based on the condition
  new_column_name <- paste0(condition_value, "_change_in_usage")
  df <- dplyr::rename(df, !!new_column_name := change_in_usage)
  
  new_column_name <- paste0(condition_value, "_prox_to_dist_shift")
  df <- dplyr::rename(df, !!new_column_name := prox_to_dist_shift) %>% dplyr::select(-condition)
  
  return(df)
})

# Optionally, name the list elements based on the file names (without the .csv extension)
names(organoid_files) <- gsub("\\_sites_with_positional_info.csv$", "", csv_files)

#make average DMAi-change in usage from lung cell lines to counter cell line specific effects
lung_DMAi_file_list = list(DMAi_files[[2]],DMAi_files[[5]])

# Function to perform inner join on chr, start, end, and strand
merge_dfs <- function(df1, df2) {
  inner_join(df1, df2, by = c("gene_name","feature_id"))
}

# Reduce the list of dataframes by performing inner joins
merged_df <- Reduce(merge_dfs, lung_DMAi_file_list)

# Calculate the average of change_in_usage columns
merged_df <- merged_df %>%
  rowwise() %>%
  mutate(avg_lung_line_change_in_usage = mean(c_across(contains("change_in_usage")), na.rm = TRUE)) 

##make average DMAi-change in usage from all tumour organoids
# Reduce the list of organoid dataframes by performing inner joins
merged_organoid_df <- Reduce(merge_dfs, organoid_files)

# Calculate the average of change_in_usage columns
merged_organoid_df <- merged_organoid_df %>%
  rowwise() %>%
  mutate(avg_organoid_change_in_usage = mean(c_across(contains("change_in_usage")), na.rm = TRUE)) 

#organoid vs lung cell line
all_APA_events_in_both_organoid_and_DMAi_average=inner_join(merged_organoid_df,merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_organoid_and_DMAi_average, x = "avg_lung_line_change_in_usage", y = "avg_organoid_change_in_usage", 
          add = "reg.line", conf.int = TRUE,  cor.coef = TRUE, 
          cor.method = "pearson", 
          xlab = "Change in PAS usage: DMAi lung tumour organoids avg", ylab = "Change in PAS usage: DMAi lung cell lines avg",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("DMAi APA - Lung tumour organoid vs Lung cancer line"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 18),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.6,0.6))+
  xlim(c(-0.56,0.57))


# Create the base plot with hexagonal binning
ggplot(all_APA_events_in_both_organoid_and_DMAi_average, aes(x = avg_lung_line_change_in_usage, y = avg_organoid_change_in_usage)) +
  geom_hex(bins = 90, aes(fill = ..count..)) +  # Use hex bins to represent density
  scale_fill_viridis_c(option = "C") +  # Apply the viridis color scale
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the regression line
  labs(
    x = "Change in pA site usage:\nDMAi lung tumour organoid average",
    y = "Change in pA site usage:\nDMAi lung cell line average"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.position = "right",legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  ylim(c(-0.6, 0.6)) +
  xlim(c(-0.6, 0.6)) +
  annotate("text", x = -0.35, y = 0.45, label = "r = 0.75", fontface = "bold",size = 10, color = "red")

###murine T cell analysis
#setwd for murine T cells
setwd("Data/Figure_2_S5/APA_CSV_files/murine_T_cells/positional_info/")
csv_files <- list.files(pattern = "top2*")  # This regex matches files that end with .csv

# Read each CSV file and store them in a list called 'DMAi_files'
murine_T_files <- lapply(csv_files, function(file) {
  df <- read.csv(file)  # Or use read_csv(file) for faster loading with readr
  
  # Select necessary columns
  df <- df %>%
    dplyr::select(gene_name, feature_id, change_in_usage, condition, prox_to_dist_shift, twostep_feature_padj)
  
  # Assume 'condition' is the same across all rows, take the first one
  condition_value <- unique(df$condition)[1]
  
  # Rename 'change_in_usage' based on the condition
  new_column_name <- paste0(condition_value, "_change_in_usage")
  df <- dplyr::rename(df, !!new_column_name := change_in_usage)
  
  new_column_name <- paste0(condition_value, "_prox_to_dist_shift")
  df <- dplyr::rename(df, !!new_column_name := prox_to_dist_shift)
  
  new_column_name <- paste0(condition_value, "_twostep_feature_padj")
  df <- dplyr::rename(df, !!new_column_name := twostep_feature_padj) %>% dplyr::select(-condition)
  
  return(df)
})

# Optionally, name the list elements based on the file names (without the .csv extension)
names(murine_T_files) <- gsub("\\_sites_with_positional_info.csv$", "", csv_files)

#Naive vs Polarised against Polarised vs Polarised + DMAi
all_APA_events_in_both_N_vs_P_and_P_vs_Pdmai=inner_join(murine_T_files[[2]],murine_T_files[[3]],by=c('gene_name','feature_id'))

#get inverse of Polarised vs Polarised + DMAi change in usage
all_APA_events_in_both_N_vs_P_and_P_vs_Pdmai$Polarised_DMSO_vs_DMAi_change_in_usage_inverted = all_APA_events_in_both_N_vs_P_and_P_vs_Pdmai$Polarised_DMSO_vs_DMAi_change_in_usage * -1

ggscatter(all_APA_events_in_both_N_vs_P_and_P_vs_Pdmai, x = "Naive_vs_Polarised_change_in_usage", y = "Polarised_DMSO_vs_DMAi_change_in_usage_inverted", 
          add = "reg.line", conf.int = TRUE,  cor.coef = TRUE, 
          cor.method = "pearson", 
          xlab = "Change in PAS usage: Naive vs Polarised", ylab = "Change in PAS usage: Polarised+DMAi vs Polarised",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("DMAi APA - Naive to Polarised vs Polarised+DMAi to Polarised"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 18),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.6,0.6))+
  xlim(c(-0.56,0.57))

# Create the base plot with hexagonal binning
ggplot(all_APA_events_in_both_N_vs_P_and_P_vs_Pdmai, aes(x = Naive_vs_Polarised_change_in_usage, y = Polarised_DMSO_vs_DMAi_change_in_usage_inverted)) +
  geom_hex(bins = 90, aes(fill = ..count..)) +  # Use hex bins to represent density
  scale_fill_viridis_c(option = "C") +  # Apply the viridis color scale
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the regression line
  labs(
    x = "Change in PAS usage: Naive vs Activated",
    y = "Change in PAS usage:\nActivated DMAi vs Activated"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 16),
    legend.position = "right",legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  ylim(c(-0.52, 0.52)) +
  xlim(c(-0.52, 0.52)) +
  annotate("text", x = -0.35, y = 0.45, label = "r = 0.5", fontface = "bold",size = 10, color = "red")
