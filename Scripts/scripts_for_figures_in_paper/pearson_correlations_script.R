library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)

# Set the working directory to the directory containing the CSV files
# (Adjust the path according to your specific directory structure)
setwd("/path/to/your/directory")
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/panel_of_cancer_lines/positional_info")

# List top2 files in the directory
csv_files <- list.files(directory, pattern = "*top2*")  # This regex matches files that have top2 in their name

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

#pearson plot: MCF7 vs LU99
all_APA_events_in_both_LU99_and_MCF7=inner_join(DMAi_files[[4]],DMAi_files[[5]],by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_LU99_and_MCF7, x = "LU99_change_in_usage", y = "MCF7_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson", 
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots
          xlab = "change in PAS usage: LU99", ylab = "change in PAS usage: MCF7")+
  labs(title = paste("LU99 vs MCF7"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  ylim(c(-0.8,0.8))+
  xlim(c(-0.8,0.8))+
  theme_bw() +  # Uses a minimal theme for clean aesthetics
  theme(
    plot.title = element_text(hjust = 0.5, size = 34, face = "bold"),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 28),
    legend.position = "none"  # Removes legend since there are no groups
  ) 

# Create the base plot with hexagonal binning
ggplot(all_APA_events_in_both_LU99_and_MCF7, aes(x = LU99_change_in_usage, y = MCF7_change_in_usage)) +
  geom_hex(bins = 90, aes(fill = ..count..)) +  # Use hex bins to represent density
  scale_fill_viridis_c(option = "C") +  # Apply the viridis color scale
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the regression line
  labs(
    x = "Change in PAS usage: MCF7",
    y = "Change in PAS usage: LU99",
    title = "DMAi APA - MCF7 vs LU99"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 22),
    legend.position = "right"
  ) +
  ylim(c(-0.6, 0.6)) +
  xlim(c(-0.6, 0.6)) +
  annotate("text", x = -0.4, y = 0.45, label = "r = 0.71", fontface = 'bold', size = 13, color = "red")

#pearson plot: PANC0403 vs GP2D
all_APA_events_in_both_PANC0403_and_GP2D=inner_join(DMAi_files[[6]],DMAi_files[[1]],by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_PANC0403_and_GP2D, x = "PANC0403_change_in_usage", y = "GP2D_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson", 
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots
          xlab = "change in PAS usage: PANC0403", ylab = "change in PAS usage: GP2D")+
  labs(title = paste("PANC0403 vs GP2D"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  ylim(c(-0.8,0.8))+
  xlim(c(-0.8,0.8))+
  theme_bw() +  # Uses a minimal theme for clean aesthetics
  theme(
    plot.title = element_text(hjust = 0.5, size = 34, face = "bold"),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 28),
    legend.position = "none"  # Removes legend since there are no groups
  )

# Create the base plot with hexagonal binning
ggplot(all_APA_events_in_both_PANC0403_and_GP2D, aes(x = PANC0403_change_in_usage, y = GP2D_change_in_usage)) +
  geom_hex(bins = 90, aes(fill = ..count..)) +  # Use hex bins to represent density
  scale_fill_viridis_c(option = "C") +  # Apply the viridis color scale
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the regression line
  labs(
    x = "Change in PAS usage: PANC0403",
    y = "Change in PAS usage: GP2D",
    title = "DMAi APA - GP2D vs PANC0403"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 22),
    legend.position = "right"
  ) +
  ylim(c(-0.6, 0.6)) +
  xlim(c(-0.6, 0.6)) +
  annotate("text", x = -0.4, y = 0.45, label = "r = 0.66", fontface = 'bold', size = 13, color = "red")

#pearson plot: H838 vs SUM149PT
all_APA_events_in_both_H838_and_SUM149PT=inner_join(DMAi_files[[2]],DMAi_files[[7]],by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_H838_and_SUM149PT, x = "H838_change_in_usage", y = "SUM149PT_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson", 
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots
          xlab = "change in PAS usage: H838", ylab = "change in PAS usage: SUM149PT")+
  labs(title = paste("H838 vs SUM149PT"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  ylim(c(-0.8,0.8))+
  xlim(c(-0.8,0.8))+
  theme_bw() +  # Uses a minimal theme for clean aesthetics
  theme(
    plot.title = element_text(hjust = 0.5, size = 34, face = "bold"),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 28),
    legend.position = "none"  # Removes legend since there are no groups
  ) 

# Create the base plot with hexagonal binning
ggplot(all_APA_events_in_both_H838_and_SUM149PT, aes(x = H838_change_in_usage, y = SUM149PT_change_in_usage)) +
  geom_hex(bins = 90, aes(fill = ..count..)) +  # Use hex bins to represent density
  scale_fill_viridis_c(option = "C") +  # Apply the viridis color scale
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the regression line
  labs(
    x = "Change in PAS usage: H838",
    y = "Change in PAS usage: SUM149PT",
    title = "DMAi APA - H838 vs SUM149PT"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 22),
    legend.position = "right"
  ) +
  ylim(c(-0.6, 0.6)) +
  xlim(c(-0.6, 0.6)) +
  annotate("text", x = -0.4, y = 0.45, label = "r = 0.62", fontface = 'bold', size = 13, color = "red")

#pearson plot: HCT116 p53 +/+ vs U2OS
all_APA_events_in_both_HCT116_p53plus_and_U2OS=inner_join(DMAi_files[[3]],DMAi_files[[8]],by=c('gene_name','feature_id'))
ggscatter(all_APA_events_in_both_HCT116_p53plus_and_U2OS, x = "HCT116_p53plus_change_in_usage", y = "U2OS_change_in_usage", 
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
          cor.method = "pearson", 
          shape = 19,            # Specifies the shape of points; 19 is solid circle
          color = "#2B496A",    # Specifies the color of points
          dot.size = 3,          # Size of the dots
          xlab = "change in PAS usage: HCT116 p53 +/+", ylab = "change in PAS usage: U2OS")+
  labs(title = paste("HCT116 p53 +/+ vs U2OS"))+
  theme(plot.title = element_text(hjust = 0.5))  +
  ylim(c(-0.8,0.8))+
  xlim(c(-0.8,0.8))+
  theme_bw() +  # Uses a minimal theme for clean aesthetics
  theme(
    plot.title = element_text(hjust = 0.5, size = 34, face = "bold"),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 28),
    legend.position = "none"  # Removes legend since there are no groups
  ) 

# Create the base plot with hexagonal binning
ggplot(all_APA_events_in_both_HCT116_p53plus_and_U2OS, aes(x = HCT116_p53plus_change_in_usage, y = U2OS_change_in_usage)) +
  geom_hex(bins = 90, aes(fill = ..count..)) +  # Use hex bins to represent density
  scale_fill_viridis_c(option = "C") +  # Apply the viridis color scale
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add the regression line
  labs(
    x = "Change in PAS usage: HCT116 p53 +/+",
    y = "Change in PAS usage: U2OS",
    title = "DMAi APA - HCT116 p53 +/+ vs U2OS"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 22),
    legend.position = "right"
  ) +
  ylim(c(-0.6, 0.6)) +
  xlim(c(-0.6, 0.6)) +
  annotate("text", x = -0.4, y = 0.45, label = "r = 0.65", fontface = 'bold', size = 13, color = "red")

### lung cell line vs organoid comparison
#setwd for organoids
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/patient_derived_organoids/positional_info")

# List top2 files in the directory
csv_files <- list.files(directory, pattern = "*top2*")  # This regex matches files that have top2 in their name

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
lung_DMAi_file_list = list(DMAi_files[[2]],DMAi_files[[4]])

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


### DMAi vs non-DMAi datasets analysis
#setwd for non DMAi
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig4/non_DMAi_datasets")

# List top2 files in the directory
csv_files <- list.files(directory, pattern = "*top2*")  # This regex matches files that have top2 in their name

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

#siCFIM25 vs DMAi (calculated r = -0.51)
all_APA_events_in_both_siCFIM25_and_DMAi_average=inner_join(non_DMAi_files[[7]],merged_df,by=c('gene_name','feature_id'))
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
all_APA_events_in_both_siELAVL1_and_DMAi_average=inner_join(non_DMAi_files[[8]],merged_df,by=c('gene_name','feature_id'))
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
all_APA_events_in_both_siPCF11_and_DMAi_average=inner_join(non_DMAi_files[[9]],merged_df,by=c('gene_name','feature_id'))
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
all_APA_events_in_both_siSam68_and_DMAi_average=inner_join(non_DMAi_files[[10]],merged_df,by=c('gene_name','feature_id'))
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
all_APA_events_in_both_siSUPT4H1_and_DMAi_average=inner_join(non_DMAi_files[[11]],merged_df,by=c('gene_name','feature_id'))
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
all_APA_events_in_both_siXRN2_and_DMAi_average=inner_join(non_DMAi_files[[12]],merged_df,by=c('gene_name','feature_id'))
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

###murine T cell analysis
#setwd for murine T cells
setwd("../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig2/murine_T_cells")

# List top2 files in the directory
csv_files <- list.files(directory, pattern = "*top2*")  # This regex matches files that have top2 in their name

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
all_APA_events_in_both_N_vs_P_and_P_vs_Pdmai=inner_join(murine_T_files[[1]],murine_T_files[[2]],by=c('gene_name','feature_id'))

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
