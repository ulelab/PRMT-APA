#load libraries
library(tidyverse)
library(ggplot2)

#setwd
setwd('/Users/llywelyngriffith/Documents/AZ_postdoc/3_end_seq_timecourse_LU99/Nobby_APA_analysis/common_atlas/dedup/github/metascape/drimseq_filter_background')

#read in file
sig_24hrs = read.csv('sig_24hrs_metascape_with_drimseq_background/sig_24hr_result_for_bubble_plot.csv')

#unite category column
sig_24hrs = sig_24hrs %>%  unite("GO Category", GO, Description, sep = ":")

#create -log10(p) column
sig_24hrs = sig_24hrs %>%
  mutate(`-log10(p)` = -1 * Log10.P.) %>%
  select(-Log10.P.)


# Calculate the maximum -log10(p) per GO Category and arrange it
sig_24hrs <- sig_24hrs %>%
  group_by(`GO Category`) %>%
  mutate(max_log10p = max(`-log10(p)`)) %>%
  ungroup() %>%
  arrange(desc(max_log10p)) %>%
  mutate(`GO Category` = fct_reorder(`GO Category`, max_log10p, .desc = FALSE))

#make bubble plot
# Custom function to split long labels
split_label <- function(x) {
  sapply(strwrap(x, width = 50, simplify = FALSE), paste, collapse = "\n")
}

ggplot(sig_24hrs, aes(x = `-log10(p)`, y = `GO Category`, size = Count, color = `-log10(p)`)) +
  geom_point(alpha = 1, show.legend = TRUE) +
  scale_color_gradient(low = "#face45", high = "#c80900") +
  labs(
    x = "-log10(p)",
    y = "GO Category",
    size = "Count",
    color = "-log10(p) Value"
  ) +
  scale_size(range = c(5, 15)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 17, angle = 45, hjust = 1),
    axis.title = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  scale_y_discrete(labels = split_label)

