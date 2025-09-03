library(tidyverse)
library(ggplot2)
library(ggpattern)

# Specify the directory
results_dir <- "/Users/k2362866/Documents/AZ_postdoc/files_for_charlotte/github/motif_and_clip_analysis/CFIM_mitigation_categories/SRA_coverage_results/shCtrl_K562_HepG2_CFIM_mitigation_categories/"

# List all files in the directory that end with 'results.csv'
results_files <- list.files(
  path = results_dir, 
  pattern = "results\\.csv$", 
  full.names = TRUE
)

# Read each file into a list
results_list <- lapply(results_files, function(file) {
  read.csv(file, header = TRUE, stringsAsFactors = FALSE)
})

# Combine all data frames in results_list into one table
combined_results <- do.call(rbind, results_list)

#refactor
combined_results$condition <- factor(
  combined_results$condition,
  levels  = c("DMAi_no_mit", "DMAi_part_mit", "DMAi_full_mit"),          # original values
  labels  = c("DMAi no mitigation",
              "DMAi part mitigation",
              "DMAi full mitigation")                     # what you want to see
)

combined_results$condition <- factor(
  combined_results$condition,
  levels = c("DMAi no mitigation",
             "DMAi part mitigation",
             "DMAi full mitigation")
)

#split into HepG2 and K562
HepG2_results = combined_results %>% filter(sample == 'shCtrl-HepG2-R1' | sample == 'shCtrl-HepG2-R2') %>% filter(condition != 'JTE607' & condition != 'JTE607_ctrl' & condition != 'siCFIM25')
K562_results = combined_results %>% filter(sample == 'shCtrl-K562-R1' | sample == 'shCtrl-K562-R2')  %>% filter(condition != 'JTE607' & condition != 'JTE607_ctrl' & condition != 'siCFIM25')

# Create the boxplot
ggplot(combined_results,
       aes(x = condition, y = proximal_use, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(                    # <- change here
    values = c(
      "DMAi full mitigation" = "#ff594e",
      "DMAi part mitigation" = "#ff9e4e",
      "DMAi no mitigation"   = "#469de1"
    )
  ) +
  labs(
    title = "DMAi mitigation categories – HepG2 RNA-seq pPA usage estimate",
    x     = "DMAi mitigation category",
    y     = "Proximal pA usage"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title   = element_text(size = 14)
  )

ggplot(K562_results,
       aes(x = condition, y = proximal_use, fill =sample)) +
  geom_boxplot() +
  labs(
    title = "DMAi mitigation categories – shCtrl K562 RNA-seq pPA usage estimate",
    x     = "Transcript Group",
    y     = "Baseline Proximal pA usage"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title   = element_text(size = 14)
  )

ggplot(
  HepG2_results,
  aes(x = condition, y = proximal_use,
      fill = condition,            # fill by condition (your 3 groups)
      pattern = sample,            # pattern by sample (R1 vs R2)
      group = interaction(condition, sample))
) +
  geom_boxplot_pattern(
    position = position_dodge2(preserve = "single", padding = 0.2),
    outlier.shape      = NA,
    pattern_fill       = NA,       # let the fill colour show
    pattern_alpha      = 0.9,      # lighter pattern
    pattern_density    = 0.2,      # fewer stripes/hatches
    pattern_spacing    = 0.06,     # more space between lines
    pattern_colour     = "black", # softer lines
    pattern_size       = 0.2,      # thinner lines
    pattern_key_scale_factor = 0.6 # smaller legend pattern
  ) +
  scale_pattern_manual(values = c(
    "shCtrl-HepG2-R1" = "stripe",
    "shCtrl-HepG2-R2" = "crosshatch"
  )) +
  scale_fill_manual(values = c(
    "DMAi full mitigation" = "#ff594e",
    "DMAi part mitigation" = "#ff9e4e",
    "DMAi no mitigation"   = "#469de1"
  )) +
  labs(
    title   = "DMAi mitigation categories – shCtrl K562 RNA-seq pPA usage estimate",
    x       = "Transcript Group",
    y       = "Baseline Proximal pA usage",
    fill    = "Category",
    pattern = "RNA-seq Sample"
  ) +
  guides(
    fill    = guide_legend(order = 1),
    pattern = guide_legend(order = 2)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title   = element_text(size = 14, hjust = 0.2)
  )
