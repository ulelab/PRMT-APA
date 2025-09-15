library(rtracklayer)
library(GenomicRanges)
library(Matrix)
library(genomation)
library(tidyverse)
library(data.table)
library(zoo)
library(RColorBrewer)
library(ggthemes)
library(ggrepel)


#import bigwigs, for this script they need to end with forward.bigWig or reverse.bigWig
rnaseq_bw = list.files("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/Charlotte_metaprofiles/RNA_seq_files/ctrl/",full.names = TRUE, pattern = "\\.bigWig$")
experiment_name = "shCtrl_K562"

granges_list <- lapply(rnaseq_bw, function(file_path) {
  # Determine the strand from the file name
  strand <- ifelse(grepl("pos", file_path), "+", "-")
  # Import the BigWig file
  gr_data <- rtracklayer::import(file_path, format = "BigWig")
  # Assign the strand
  strand(gr_data) <- strand
  gr_data
})

# Updated merge function to maintain scores
merge_strands <- function(pos, neg) {
  if (is.null(pos) || is.null(neg)) {
    return(ifelse(is.null(pos), neg, pos))
  }
  # Combine forward and reverse
  combined <- c(pos, neg)
  return(combined)
}

sample_names <- gsub("^(.*)(pos|neg)\\.bigWig$", "\\1", basename(rnaseq_bw))
sample_names <- gsub("(pos|neg)$", "", sample_names)  # Ensure to remove the strand part

# Group GRanges objects by sample
sample_groups <- split(granges_list, sample_names)

# Process in smaller batches
batch_process <- function(batch_data) {
  lapply(batch_data, function(gr_list) {
    merge_strands(gr_list[[1]], gr_list[[2]])
  })
}

# divide your sample_groups into manageable chunks
batches <- split(sample_groups, ceiling(seq_along(sample_groups) / 2))  # Adjust the denominator for batch size
processed_batches <- lapply(batches, batch_process)

# Combine processed data
merged_samples <- do.call(c, processed_batches)
length(merged_samples)

# Normalise by library size
# Function to normalize scores by library size (CPM)
normalize_scores <- function(gr) {
  total_score <- sum(mcols(gr)$score)
  cpm <- (mcols(gr)$score / total_score) * 1e6  # Convert to counts per million
  mcols(gr)$score <- cpm  # Assign normalized scores back to the GRanges object
  return(gr)
}

# Apply normalization to each GRanges object in the list
normalized_granges_list <- lapply(merged_samples, normalize_scores)

#### read in poly A sites
#read in files
DMAi_down_proximal_sites = fread("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/bed/test_triple_removal/expression_matched/2000/common_DMAi/common_DMAi_sig_down_TUTR_proximal_sites.bed") %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    strand = V6) %>%
  unique()


non_DMAi_down_proximal_sites = fread("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/bed/test_triple_removal/expression_matched/2000/non_DMAi/non_DMAi_sig_down_TUTR_proximal_sites.bed") %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    strand = V6) %>%
  unique()

ctrl_proximal_sites = fread("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/bed/test_triple_removal/expression_matched/2000/common_DMAi/control_TUTR_proximal_sites.bed") %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    strand = V6) %>%
  unique()

DMAi_up_distal_sites = fread("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/bed/test_triple_removal/expression_matched/2000/common_DMAi/common_DMAi_sig_up_TUTR_distal_sites.bed") %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    strand = V6) %>%
  unique()


non_DMAi_up_distal_sites = fread("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/bed/test_triple_removal/expression_matched/2000/non_DMAi/non_DMAi_sig_up_TUTR_distal_sites.bed") %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    strand = V6) %>%
  unique()

ctrl_distal_sites = fread("/Users/k2362866/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/github/bed/test_triple_removal/expression_matched/2000/common_DMAi/control_TUTR_distal_sites.bed") %>%
  dplyr::rename(
    chr = V1, 
    start = V2,
    end = V3, 
    strand = V6) %>%
  unique()


nrow(DMAi_up_distal_sites)
nrow(non_DMAi_up_distal_sites)
nrow(ctrl_proximal_sites)
nrow(DMAi_down_proximal_sites)
nrow(non_DMAi_down_proximal_sites)
nrow(ctrl_distal_sites)

DMAi_down_prox = makeGRangesFromDataFrame(DMAi_down_proximal_sites, keep.extra.columns = TRUE)
non_DMAi_down_prox = makeGRangesFromDataFrame(non_DMAi_down_proximal_sites, keep.extra.columns = TRUE)
DMAi_up_dist = makeGRangesFromDataFrame(DMAi_up_distal_sites, keep.extra.columns = TRUE)
non_DMAi_up_dist = makeGRangesFromDataFrame(non_DMAi_up_distal_sites, keep.extra.columns = TRUE)
ctrl_prox = makeGRangesFromDataFrame(ctrl_proximal_sites, keep.extra.columns = TRUE)
ctrl_dist = makeGRangesFromDataFrame(ctrl_distal_sites, keep.extra.columns = TRUE)

DMAi_down_prox_2000 = resize(DMAi_down_prox, width = 2000, fix = "center")
non_DMAi_down_prox_2000 = resize(non_DMAi_down_prox, width = 2000, fix = "center")
DMAi_up_dist_2000 = resize(DMAi_up_dist, width = 2000, fix = "center")
non_DMAi_up_dist_2000 = resize(non_DMAi_up_dist, width = 2000, fix = "center")
ctrl_prox_2000 = resize(ctrl_prox, width = 2000, fix = "center")
ctrl_dist_2000 = resize(ctrl_dist, width = 2000, fix = "center")

DMAi_prox = genomation::ScoreMatrixList(target = normalized_granges_list, windows = DMAi_down_prox_2000, strand.aware=TRUE, weight.col= "score", cores=4)
DMAi_dist = genomation::ScoreMatrixList(target = normalized_granges_list, windows = DMAi_up_dist_2000,strand.aware=TRUE, weight.col= "score", cores=4 )
non_DMAi_prox = genomation::ScoreMatrixList(target = normalized_granges_list, windows = non_DMAi_down_prox_2000,strand.aware=TRUE, weight.col= "score", cores=4)
non_DMAi_dist = genomation::ScoreMatrixList(target = normalized_granges_list, windows = non_DMAi_up_dist_2000,strand.aware=TRUE, weight.col= "score", cores=4)
ctrl_prox = genomation::ScoreMatrixList(target = normalized_granges_list, windows = ctrl_prox_2000,strand.aware=TRUE, weight.col= "score", cores=4)
ctrl_dist = genomation::ScoreMatrixList(target = normalized_granges_list, windows = ctrl_dist_2000,strand.aware=TRUE, weight.col= "score", cores=4)

### ggplot
make_meta_dataframe = function(df, condition, location, sample){
  new_df = df %>%
    pivot_longer(cols=everything()) %>%
    group_by(name) %>%
    mutate(mean=mean(value), 
           name=as.numeric(gsub("V","",name)),
           condition=condition,
           location=location,
           sample=sample) %>%
    dplyr::select(name, mean, condition, location, sample) %>%
    unique()
  return(new_df)
}

combine_dataframes = function(scorematrixlist, condition, location){
  df_list = list()
  for (i in 1:length(scorematrixlist)){
    df_list[[i]] = make_meta_dataframe(as.data.frame(scorematrixlist[[i]]), condition, location,unique(sample_names)[[i]])
  }
  for_plot=bind_rows(df_list)
  return(for_plot)
}

DMAi_prox_df = combine_dataframes(DMAi_prox,"DMAi","proximal")
DMAi_dist_df = combine_dataframes(DMAi_dist,"DMAi","distal")
non_DMAi_prox_df = combine_dataframes(non_DMAi_prox,"non_DMAi","proximal")
non_DMAi_dist_df = combine_dataframes(non_DMAi_dist,"non_DMAi","distal")
ctrl_prox_df = combine_dataframes(ctrl_prox,"control","proximal")
ctrl_dist_df = combine_dataframes(ctrl_dist,"control","distal")

final_to_plot = bind_rows(DMAi_prox_df, DMAi_dist_df,non_DMAi_prox_df,non_DMAi_dist_df,ctrl_prox_df,ctrl_dist_df) 

final_to_plot = final_to_plot %>%
  mutate(location=factor(location, levels=c("proximal","distal")),
         condition=factor(condition, levels=c("DMAi","non_DMAi","control")))

ggplot(final_to_plot,aes(x=name,y=mean,color=condition)) +
  geom_line(aes(y=rollmean(mean, 7, na.pad=TRUE))) + 
  facet_wrap(~sample + location) + theme_bw()  +
  scale_color_manual(
    values = c(
      "DMAi" = "#e53b3b",
      "non_DMAi" = "#669b5f",
      "ctrl" = "darkgrey"
    ),
    labels = c(
      "DMAi",
      "non_DMAi",
      "control"
    )
  )+
  geom_vline(xintercept=c(1000), linetype="dotted") +
  scale_x_continuous(breaks = seq(0, 2000, by = 200),   # X-axis breaks every 100 units
                     labels = seq(-1000, 1000, by = 200),  # X-axis labels every 100 units
                     limits = c(0, 2000)) +
  ylab("Mean CPM") +
  xlab("Position") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 24),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 20)
    )
  

