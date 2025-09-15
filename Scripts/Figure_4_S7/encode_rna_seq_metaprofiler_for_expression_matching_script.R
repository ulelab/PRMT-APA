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
# Read in files and create GRanges objects
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

# List of file paths
bed_file_paths <- list(
  full_mit_proximal = "/Users/k2362866/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/expression_matched/2000/DMAi_full_mit_proximal_sites.bed",
  part_mit_proximal = "/Users/k2362866/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/expression_matched/2000/DMAi_part_mit_proximal_sites.bed",
  no_mit_proximal   = "/Users/k2362866/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/expression_matched/2000/no_mit_proximal_sites.bed",
  full_mit_distal   = "/Users/k2362866/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/expression_matched/2000/DMAi_full_mit_distal_sites.bed",
  part_mit_distal   = "/Users/k2362866/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/expression_matched/2000/DMAi_part_mit_distal_sites.bed",
  no_mit_distal     = "/Users/k2362866/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/common_atlas/dedup/github/siCFIM25_DMSO_vs_DMAi/bed/expression_matched/2000/no_mit_distal_sites.bed"
)

# Read and prepare GRanges objects
gr_list <- lapply(bed_file_paths, read_and_prepare_bed)

# Resize GRanges objects
resize_gr <- function(gr) {
  resize(gr, width = 2000, fix = "center")
}

gr_resized_list <- lapply(gr_list, resize_gr)

full_mit_prox = genomation::ScoreMatrixList(target = normalized_granges_list, windows = gr_resized_list[[1]], strand.aware=TRUE, weight.col= "score", cores=4)
part_mit_prox = genomation::ScoreMatrixList(target = normalized_granges_list, windows = gr_resized_list[[2]], strand.aware=TRUE, weight.col= "score", cores=4 )
no_mit_prox = genomation::ScoreMatrixList(target = normalized_granges_list, windows = gr_resized_list[[3]], strand.aware=TRUE, weight.col= "score", cores=4)
full_mit_dist = genomation::ScoreMatrixList(target = normalized_granges_list, windows = gr_resized_list[[4]], strand.aware=TRUE, weight.col= "score", cores=4)
part_mit_dist = genomation::ScoreMatrixList(target = normalized_granges_list, windows = gr_resized_list[[5]], strand.aware=TRUE, weight.col= "score", cores=4 )
no_mit_dist = genomation::ScoreMatrixList(target = normalized_granges_list, windows = gr_resized_list[[6]], strand.aware=TRUE, weight.col= "score", cores=4)

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

full_mit_prox_df = combine_dataframes(full_mit_prox,"DMAi Full Mitigation","proximal")
part_mit_prox_df = combine_dataframes(part_mit_prox,"DMAi Part Mitigation","proximal")
no_mit_prox_df = combine_dataframes(no_mit_prox,"DMAi No Mitigation","proximal")

full_mit_dist_df = combine_dataframes(full_mit_dist,"DMAi Full Mitigation","distal")
part_mit_dist_df = combine_dataframes(part_mit_dist,"DMAi Part Mitigation","distal")
no_mit_dist_df = combine_dataframes(no_mit_dist,"DMAi No Mitigation","distal")


final_to_plot = bind_rows(full_mit_prox_df, part_mit_prox_df, no_mit_prox_df,
                          full_mit_dist_df, part_mit_dist_df, no_mit_dist_df) 

final_to_plot = final_to_plot %>%
  mutate(location=factor(location, levels=c("proximal","distal")),
         condition=factor(condition, levels=c('DMAi No Mitigation','DMAi Part Mitigation','DMAi Full Mitigation')))

ggplot(final_to_plot,aes(x=name,y=mean,color=condition)) +
  geom_line(aes(y=rollmean(mean, 7, na.pad=TRUE))) + 
  facet_wrap(~sample + location) + theme_bw()  +
  scale_color_manual(
    values = c(
      "DMAi No Mitigation" = "#469de1",
      "DMAi Part Mitigation" = "#ff9e4e",
      "DMAi Full Mitigation" = "#ff594e"
    ),
    labels = c(
      "DMAi No Mitigation",
      "DMAi Part Mitigation",
      "DMAi Full Mitigation"
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
