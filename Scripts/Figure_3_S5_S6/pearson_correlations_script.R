library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)

# Set the working directory to the directory containing the CSV files
# (Adjust the path according to your specific directory structure)
setwd("/path/to/your/directory")
setwd("../../Data/Figure_2_S4/APA_CSV_files/panel_of_cancer_lines/positional_info")

# List all CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")  # This regex matches files that end with .csv

# Read each CSV file and store them in a list called 'DMAi_files'
DMAi_files <- lapply(csv_files, function(file) {
  df <- read.csv(file)  # Or use read_csv(file) for faster loading with readr
  
  # Select necessary columns
  df <- df %>%
    dplyr::select(gene_name, feature_id, change_in_usage, condition, prox_to_dist_shift)
  
  # Assume 'condition' is the same across all rows, take the first one
  condition_value <- unique(df$condition)[1]
  
  # Rename 'change_in_usage' and 'prox_to_dist_shift' based on the condition
  new_column_name <- paste0(condition_value, "_change_in_usage")
  df <- dplyr::rename(df, !!new_column_name := change_in_usage)
  
  new_column_name <- paste0(condition_value, "_prox_to_dist_shift")
  df <- dplyr::rename(df, !!new_column_name := prox_to_dist_shift) %>% dplyr::select(-condition)
  
  return(df)
})

# Optionally, name the list elements based on the file names (without the .csv extension)
names(DMAi_files) <- gsub("\\_sites_with_positional_info.csv$", "", csv_files)

### DMAi vs non-DMAi datasets analysis
#setwd for non DMAi
setwd("../../Data/Figure_3_S5_S6/APA_CSV_files/non_DMAi_datasets/positional_info/")
csv_files <- list.files(pattern = "\\.csv$")  # This regex matches files that end with .csv

# Read each CSV file and store them in a list called 'non_DMAi_files'
non_DMAi_files <- lapply(csv_files, function(file) {
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
names(non_DMAi_files) <- gsub("\\_sites_with_positional_info.csv$", "", csv_files)

#make average DMAi-change in usage from 4 most sensitive cell lines to counter cell line specific effects
DMAi_file_list = list(DMAi_files[[1]],DMAi_files[[3]],DMAi_files[[4]],DMAi_files[[5]])

# Reduce the list of dataframes by performing inner joins
merged_df <- Reduce(merge_dfs, DMAi_file_list)

# Calculate the average of change_in_usage columns
merged_df <- merged_df %>%
  rowwise() %>%
  mutate(avg_DMAi_change_in_usage = mean(c_across(contains("change_in_usage")), na.rm = TRUE)) 

#CBPi vs DMAi (calculated r = 0.029)
all_APA_events_in_both_CBPi_HepG2_and_DMAi_average=inner_join(non_DMAi_files[[1]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_CBPi_HepG2_and_DMAi_average, x = "CBPi_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - CBPi", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - CBPi vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#CDK11i vs DMAi (calculated r = 0.12)
all_APA_events_in_both_CDK11i_HepG2_and_DMAi_average=inner_join(non_DMAi_files[[2]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_CDK11i_HepG2_and_DMAi_average, x = "CDK11i_24hrs_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - CDK11i", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - CDK11i vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#CDK12i vs DMAi (calculated r = -0.069)
all_APA_events_in_both_CDK12i_HepG2_and_DMAi_average=inner_join(non_DMAi_files[[3]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_CDK12i_HepG2_and_DMAi_average, x = "CDK12i_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - CDK12i", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - CDK12i vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#CLKi vs DMAi (calculated r = 0.46)
all_APA_events_in_both_CLKi_HepG2_and_DMAi_average=inner_join(non_DMAi_files[[4]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_CLKi_HepG2_and_DMAi_average, x = "CLKi_N_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - CLKi", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - CLKi vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))


#JTE607 vs DMAi (calculated r = 0.47)
all_APA_events_in_both_JTE_HepG2_and_DMAi_average=inner_join(non_DMAi_files[[5]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_JTE_HepG2_and_DMAi_average, x = "JTE607_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - JTE607 HepG2", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - JTE607 vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

# Create the base plot with hexagonal binning
ggplot(all_APA_events_in_both_JTE_HepG2_and_DMAi_average, aes(x = JTE607_change_in_usage, y = avg_DMAi_change_in_usage)) +
  geom_hex(bins = 90, aes(fill = ..count..)) +  # Use hex bins to represent density
  scale_fill_viridis_c(option = "C") +  # Apply the viridis color scale
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the regression line
  labs(
    x = "Change in PAS usage: CPSF73i HepG2",
    y = "Change in PAS usage:\nDMAi cell line avg",
    title = "DMAi vs CPSF73i"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 20),
    legend.position = "right",legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  ylim(c(-0.65, 0.65)) +
  xlim(c(-0.65, 0.65)) +
  annotate("text", x = -0.42, y = 0.6, label = "r = 0.47", fontface = "bold",size = 10, color = "red")

#METTL3i vs DMAi (calculated r = 0.024)
all_APA_events_in_both_METTL3i_and_DMAi_average=inner_join(non_DMAi_files[[6]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_METTL3i_and_DMAi_average, x = "METTL3i_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - METTL3i", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - METTL3i vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#RBBP6 KD vs DMAi (calculated r = 0.12)
all_APA_events_in_both_RBBP6_KD_and_DMAi_average=inner_join(non_DMAi_files[[7]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_RBBP6_KD_and_DMAi_average, x = "RBBP6_kd_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - RBBP6 KD", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - RBBP6 KD vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#siCFIM25 vs DMAi (calculated r = -0.51)
all_APA_events_in_both_siCFIM25_and_DMAi_average=inner_join(non_DMAi_files[[8]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_siCFIM25_and_DMAi_average, x = "siCFIM25_DMSO_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - siCFIM25", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - siCFIM25 vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#siELAVL1 vs DMAi (calculated r = -0.2)
all_APA_events_in_both_siELAVL1_and_DMAi_average=inner_join(non_DMAi_files[[9]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_siELAVL1_and_DMAi_average, x = "siELAVL1_DMSO_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage - siELAVL1", ylab = "Change in PAS usage - DMAi cell line average",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - siELAVL1 vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#siPCF11 vs DMAi (calculated r = 0.58)
all_APA_events_in_both_siPCF11_and_DMAi_average=inner_join(non_DMAi_files[[10]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_siPCF11_and_DMAi_average, x = "siPCF11_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage: siPCF11 HeLa", ylab = "Change in PAS usage - DMAi cell line avg",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - siPCF11 vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

# Create the base plot with hexagonal binning
ggplot(all_APA_events_in_both_siPCF11_and_DMAi_average, aes(x = siPCF11_change_in_usage, y = avg_DMAi_change_in_usage)) +
  geom_hex(bins = 90, aes(fill = ..count..)) +  # Use hex bins to represent density
  scale_fill_viridis_c(option = "C") +  # Apply the viridis color scale
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the regression line
  labs(
    x = "Change in PAS usage: siPCF11 HeLa",
    y = "Change in PAS usage:\nDMAi cell line avg",
    title = "DMAi vs siPCF11"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 20),
    legend.position = "right",legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  ylim(c(-0.65, 0.65)) +
  xlim(c(-0.65, 0.65)) +
  annotate("text", x = -0.42, y = 0.6, label = "r = 0.58", fontface = "bold",size = 10, color = "red")

#siSam68 vs DMAi (calculated r = 0.21)
all_APA_events_in_both_siSam68_and_DMAi_average=inner_join(non_DMAi_files[[11]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_siSam68_and_DMAi_average, x = "siSam68_DMSO_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage: siSam68", ylab = "Change in PAS usage - DMAi cell line avg",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - siSam68 vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#siSUPT4H1 vs DMAi (calculated r = 0.07)
all_APA_events_in_both_siSUPT4H1_and_DMAi_average=inner_join(non_DMAi_files[[12]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_siSUPT4H1_and_DMAi_average, x = "siSUPT4H1_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage: siSUPT4H1", ylab = "Change in PAS usage - DMAi cell line avg",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - siSUPT4H1 vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

#siXrn2 vs DMAi (calculated r = 0.096)
all_APA_events_in_both_siXRN2_and_DMAi_average=inner_join(non_DMAi_files[[13]],merged_df,by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_siXRN2_and_DMAi_average, x = "siXRN2_change_in_usage", y = "avg_DMAi_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson",
          xlab = "Change in PAS usage: siXRN2", ylab = "Change in PAS usage - DMAi cell line avg",
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots)
          title = paste("APA Pearson r Correlation - siXRN2 vs DMAi"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"  # Removes legend since there are no groups
  ) +
  ylim(c(-0.71,0.71))+
  xlim(c(-0.71,0.71))

## make correlation plot
df_corr <- data.frame(
  dataset = c(
    "CBPi", "CDK11i", "CDK12i", "CLKi", 
    "CPSF73i", "METTL3i", "siCFIM25", "siELAVL1", 
    "siPCF11", "siKHDRBS1", "siSUPT4H1", "siXRN2"
  ),
  r_value = c(
    0.029, 0.12, -0.069, 0.46, 
    0.47, 0.024, -0.51, -0.2, 
    0.58, 0.21, 0.07, 0.096
  )
)

# 1) Order the x-axis factor levels so that highest positive correlation is on the left
#    and most negative is on the right:
df_corr$dataset <- factor(df_corr$dataset, levels = df_corr$dataset[order(df_corr$r_value, decreasing = TRUE)])

# 2) Plot a single-row "correlogram"-style chart
ggplot(df_corr, aes(x = dataset, y = "DMAi", fill = r_value, size = abs(r_value))) +
  geom_point(shape = 21, stroke = 0.8, color = "black") +  # Outline in black
  scale_size(range = c(5, 20), guide = FALSE) +
  # Use a less-white midpoint so small r-values are clearer
  scale_fill_gradient2(
    midpoint = 0,
    low = "navy",
    mid = "white",
    high = "red",
    limits = c(-0.65, 0.65)
  )+
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )  +
  labs(
    x = "APA Datasets",
    fill = "Pearson correlation (r)"
  )
