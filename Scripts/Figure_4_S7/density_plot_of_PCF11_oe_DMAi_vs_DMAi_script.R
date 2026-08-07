library(tidyverse)
library(ggplot2)
library(patchwork)

# Directory containing the CSV files
pcf11_dmai_dir <- "Data/Figure_4_S7/APA_CSV_files/APA_classified/PCF11_DMAi"

# List all CSV files
csv_files <- list.files(
  path = pcf11_dmai_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Read each CSV file and store them in a list called 'PCF11_oe_files' - need to only keep TUTR (and add prox-to-dist shift column)
PCF11_oe_files <- lapply(csv_files, function(file) {
  df <- read.csv(file)  # Or use read_csv(file) for faster loading with readr
  
  df = df %>%
    group_by(gene_name) %>%
    mutate(prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] - change_in_usage[pA_type == "proximal pA"])
    
  #dplyr::select necessary columns
  df <- df %>%
    dplyr::filter(UTR_type == 'TUTR') %>%
    dplyr::select("gene_name", "feature_id", "chr", "start", "end", "strand", "UTR_type",
                  "pA_usage_control", "pA_usage_treatment", "change_in_usage", "prox_to_dist_shift",
                  "twostep_feature_padj", "condition", "pA_type", "mean_usage")
  
  # Assume 'condition' is the same across all rows, take the first one
  condition_value <- unique(df$condition)[1]
  
  # List of columns to rename
  columns_to_rename <- c("prox_to_dist_shift", "twostep_feature_padj", "change_in_usage", 
                         "mean_usage", "pA_usage_control", "pA_usage_treatment")
  
  # Rename columns based on the condition using rename_with()
  df <- df %>%
    rename_with(.cols = all_of(columns_to_rename),
                .fn = ~ paste0(condition_value, "_", .x))
  
  return(df)
})

# Function to perform inner join on chr, start, end, and strand
merge_dfs <- function(df1, df2) {
  inner_join(df1, df2, by = c("gene_name","feature_id","chr", "start", "end", "strand","pA_type","UTR_type"))
}

# Reduce the list of dataframes by performing inner joins
merged_df <- Reduce(merge_dfs, PCF11_oe_files)

merged_df <- merged_df %>%
  group_by(gene_name) %>%
  filter(
    n() == 2,
    any(pA_type == "proximal pA"),
    any(pA_type == "distal pA")
  ) %>%
  ungroup()

#for defining dmai lengthened genes
dmai_effect_all <- function(merged_df, cl) {
  
  s <- glue::glue(
    "LU99_PCF11_{cl}_no_dox_DMSO_vs_no_dox_DMAi_prox_to_dist_shift"
  )
  
  p <- glue::glue(
    "LU99_PCF11_{cl}_no_dox_DMSO_vs_no_dox_DMAi_twostep_feature_padj"
  )
  
  merged_df %>%
    filter(pA_type == "distal pA") %>%
    transmute(
      gene_name,
      feature_id,
      shift = .data[[s]],
      padj  = .data[[p]]
    ) %>%
    dplyr::rename(
      !!paste0("dmai_shift_", cl) := shift,
      !!paste0("dmai_padj_", cl) := padj
    )
}

#set up per clone config
clone_config <- list(
  c3 = list(clone = "c3", dose = "30ng"),
  c7 = list(clone = "c7", dose = "20ng")
)

# PCF11 OE effect WITHIN DMAi (across ALL genes, not just lengthened)
pcf11_shortening_in_dmai <- function(merged_df, cl, dose) {
  s <- glue::glue("LU99_PCF11_{cl}_no_dox_DMAi_vs_{dose}_dox_DMAi_prox_to_dist_shift")
  merged_df %>%
    dplyr::filter(pA_type == "distal pA") %>%
    dplyr::transmute(gene_name, feature_id, shift = .data[[s]]) %>%
    dplyr::rename(!!paste0("pcf11_shift_", cl) := shift)
}

# Build the full table: every distal-pA gene, its DMAi effect + PCF11-in-DMAi shortening
full <- dmai_effect_all(merged_df, "c3") %>%
  
  inner_join(
    dmai_effect_all(merged_df, "c7"),
    by = c("gene_name", "feature_id")
  ) %>%
  
  inner_join(
    pcf11_shortening_in_dmai(merged_df, "c3", "30ng"),
    by = c("gene_name", "feature_id")
  ) %>%
  
  inner_join(
    pcf11_shortening_in_dmai(merged_df, "c7", "20ng"),
    by = c("gene_name", "feature_id")
  ) %>%
  
  mutate(
    
    pcf11_shift_mean =
      rowMeans(
        cbind(pcf11_shift_c3,
              pcf11_shift_c7),
        na.rm = TRUE
      ),
    
    dmai_shift_mean =
      rowMeans(
        cbind(dmai_shift_c3,
              dmai_shift_c7),
        na.rm = TRUE
      ),
    
    is_dmai_lengthened =
      !is.na(dmai_shift_c3) &
      !is.na(dmai_shift_c7) &
      dmai_shift_c3 > 0.2 &
      dmai_shift_c7 > 0.2 &
      (
        (!is.na(dmai_padj_c3) & dmai_padj_c3 < 0.05) |
          (!is.na(dmai_padj_c7) & dmai_padj_c7 < 0.05)
      )
  )

# Compare PCF11's shortening (negative shift = shortening) between groups
full %>%
  dplyr::group_by(is_dmai_lengthened) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_pcf11_shift   = mean(pcf11_shift_mean, na.rm = TRUE),
    median_pcf11_shift = median(pcf11_shift_mean, na.rm = TRUE),
    .groups = "drop"
  )

# Is PCF11 shortening stronger (more negative) at DMAi-lengthened genes?
wilcox.test(pcf11_shift_mean ~ is_dmai_lengthened, data = full,
            alternative = "greater")
# 'greater' tests FALSE > TRUE, i.e. non-lengthened shift > lengthened shift;
# since shortening is NEGATIVE, lengthened genes being MORE shortened = smaller/more negative
# → check direction and flip alternative as needed after seeing the means

ggplot(full, aes(x = pcf11_shift_mean, fill = is_dmai_lengthened)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#d73027"),
                    labels = c("Not DMAi-lengthened", "DMAi-lengthened")) +
  labs(x = "PCF11 OE shift in DMAi background\n(negative = shortening)",
       fill = NULL,
       title = "PCF11 shortening is concentrated at DMAi-lengthened genes") +
  theme_classic(base_size = 14)

# Add baseline proximal usage in the DMAi background as a covariate
prox_usage_dmai <- function(merged_df, cl, dose) {
  u <- glue::glue("LU99_PCF11_{cl}_no_dox_DMAi_vs_{dose}_dox_DMAi_pA_usage_control")
  merged_df %>%
    dplyr::filter(pA_type == "proximal pA") %>%   # proximal site usage
    dplyr::transmute(gene_name, feature_id_prox = feature_id,
                     prox_usage = .data[[u]])
}

# Build proximal-usage tables for each clone (proximal site usage in the no-dox DMAi baseline)
prox_usage_c3 <- prox_usage_dmai(merged_df, "c3", "30ng") %>%
  dplyr::rename(prox_usage_c3 = prox_usage)
prox_usage_c7 <- prox_usage_dmai(merged_df, "c7", "20ng") %>%
  dplyr::rename(prox_usage_c7 = prox_usage)

# Each has one proximal row per gene; join to full by gene_name
full <- full %>%
  dplyr::left_join(prox_usage_c3 %>% dplyr::select(gene_name, prox_usage_c3),
                   by = "gene_name") %>%
  dplyr::left_join(prox_usage_c7 %>% dplyr::select(gene_name, prox_usage_c7),
                   by = "gene_name") %>%
  dplyr::mutate(
    prox_usage = rowMeans(cbind(prox_usage_c3, prox_usage_c7), na.rm = TRUE)
  )

# Ensure the grouping factor is labelled consistently in `full`
plot_df <- full %>%
  dplyr::filter(!is.na(pcf11_shift_mean), !is.na(prox_usage)) %>%
  dplyr::mutate(
    group = factor(is_dmai_lengthened,
                   levels = c(FALSE, TRUE),
                   labels = c("Not DMAi-lengthened", "DMAi-lengthened"))
  )

grp_cols <- c("Not DMAi-lengthened" = "grey70",
              "DMAi-lengthened"     = "#d73027")

# --- Panel 1: PCF11 OE shift density (your existing plot) ---
p_shift <- ggplot(plot_df, aes(x = pcf11_shift_mean, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_manual(values = grp_cols) +
  labs(x = "PCF11 OE shift in DMAi background\n(negative = shortening)",
       y = "density", fill = NULL,
       title = "PCF11 shortening is concentrated at DMAi-lengthened genes") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

# --- Panel 2: baseline proximal usage, boxplot per group ---
p_usage <- ggplot(plot_df, aes(x = group, y = prox_usage, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, alpha = 0.1, size = 0.4) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(x = NULL, y = "Baseline proximal usage\n(DMAi, no dox)",
       title = "Baseline headroom differs between groups") +
  theme_classic(base_size = 14)

# --- Combine ---
p_shift / p_usage + plot_layout(heights = c(2, 1.3))

# Option 1: stratify by proximal-usage bin and compare within bins
full_matched <- full %>%
  dplyr::filter(!is.na(prox_usage), !is.na(pcf11_shift_mean)) %>%
  dplyr::mutate(usage_bin = ntile(prox_usage, 5)) %>%
  dplyr::group_by(usage_bin, is_dmai_lengthened) %>%
  dplyr::summarise(
    median_shift = median(pcf11_shift_mean, na.rm = TRUE),
    lower = quantile(pcf11_shift_mean, 0.25, na.rm = TRUE),
    upper = quantile(pcf11_shift_mean, 0.75, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    is_dmai_lengthened = factor(is_dmai_lengthened,
                                levels = c(FALSE, TRUE),
                                labels = c("Not DMAi-lengthened", "DMAi-lengthened"))
  )

# if lengthened > non-lengthened shortening WITHIN each usage bin, it's specific, not headroom

ggplot(full_matched,
       aes(x = factor(usage_bin), y = median_shift,
           fill = is_dmai_lengthened)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8), width = 0.2, colour = "grey30") +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8),
            vjust = ifelse(full_matched$median_shift < 0, 1.3, -0.5),
            size = 3) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(values = c("Not DMAi-lengthened" = "grey70",
                               "DMAi-lengthened"     = "#d73027")) +
  labs(
    x = "Proximal-usage quintile (headroom, low → high)",
    y = "Median PCF11 OE shift in DMAi\n(negative = shortening)",
    fill = NULL,
    title = "PCF11 shortening at DMAi-lengthened vs other genes, matched on headroom"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

# Prepare: bin, label, colour
plot_df <- full %>%
  dplyr::filter(!is.na(pcf11_shift_mean), !is.na(prox_usage)) %>%
  dplyr::mutate(
    usage_bin = ntile(prox_usage, 3),
    usage_bin = factor(
      usage_bin,
      levels = c(3,2,1),
      labels = c(
        "High proximal usage",
        "Mid proximal usage",
        "Low proximal usage"
      )
    ),
    group = factor(is_dmai_lengthened,
                   levels = c(FALSE, TRUE),
                   labels = c("Not DMAi-lengthened", "DMAi-lengthened"))
  )

grp_cols <- c("Not DMAi-lengthened" = "grey70",
              "DMAi-lengthened"     = "#d73027")

# counts per quintile per group
n_labels <- plot_df %>%
  dplyr::group_by(usage_bin, group) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

shift_tests <- plot_df %>%
  dplyr::group_by(usage_bin) %>%
  dplyr::summarise(p = wilcox.test(pcf11_shift_mean ~ group)$p.value, .groups = "drop") %>%
  dplyr::mutate(padj = p.adjust(p, method = "BH"),
                label = ifelse(padj < 2.2e-16, "p < 2e-16", paste0("p = ", signif(padj, 2))))

p_shift <- ggplot(plot_df, aes(x = pcf11_shift_mean, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  
  geom_label(
    data = shift_tests,
    aes(x = -0.9, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.4,
    size = 3.5,
    fill = "white"
  ) +
  
  facet_wrap(~ usage_bin, nrow = 1) +
  scale_fill_manual(values = grp_cols) +
  labs(
    x = "PCF11 OE prox-to-dist shift in DMAi",
    y = "density",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

# --- Bottom row: proximal usage boxplot + jittered points, faceted by quintile ---
p_usage <- ggplot(plot_df, aes(x = group, y = prox_usage, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.55) +
  geom_text(data = n_labels, aes(x = group, y = 1.05, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ usage_bin, nrow = 1) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  coord_cartesian(ylim = c(0, 1.12), clip = "off") +   # headroom for the labels
  labs(x = NULL, y = "Baseline proximal usage\n(DMAi, no dox)",
       title = "Proximal usage matched groups") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

p_shift / p_usage + plot_layout(heights = c(1.4, 1.1))
