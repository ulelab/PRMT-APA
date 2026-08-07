#read in packages
library(tidyverse)
library(data.table)

#read in file
Shechter_table = read.csv('Data/Figure_1_S1_S2_S3_S4/methylation_sites_multi_study_list_Shechter_et_al.csv')

# visualisations
study_metadata <- tribble(
  ~reference_pattern,  ~study_label,          ~model,           ~model_type,
  "Fedoriw",           "Fedoriw et al., 2019",         "DLBCL/Pancreatic\ncancer line panel", "Cancer cell line",
  "Fong",              "Fong et al., 2019",             "NB4",           "Cancer cell line",
  "Guo",               "Guo et al., 2014",              "HCT116",              "Cancer cell line",
  "Hornbeck",          "Hornbeck et al., 2015",  "PhosphoSitePlus",          "Database",
  "Larsen",            "Larsen et al., 2016",           "HEK293",            "Immortalised line",
  "Li",                "Li et al., 2021",               "HEK293",               "Immortalised line",
  "Lim",               "Lim et al., 2020",              "Colorectal cancer patient tissue",               "Patient tumour tissue",
  "Musiani",           "Musiani et al., 2019",          "HeLa",               "Cancer cell line",
  "Wei",               "Wei et al., 2020",              "PANC-1",              "Cancer cell line",
  "This study",        "Maron et al., 2021",            "A549",              "Cancer cell line"
)

# Define CPA factors of interest
cpa_factors_of_interest <- c(
  # CPSF complex
  "CPSF1", "CPSF2", "CPSF3", "CPSF4", "WDR33", "FIP1L1",
  # CFIm complex  
  "NUDT21", "CPSF6", "CPSF7",
  # CFIIm complex
  "PCF11", "CLP1",
  # CSTF complex
  "CSTF1", "CSTF2", "CSTF3", "CSTF2T",
  # Other core CPA
  "RBBP6", "SYMPK"
)

# Define complex membership for grouping in the plot
cpa_complex_map <- tribble(
  ~gene_symbol, ~complex,
  "CPSF1",  "CPSF", "CPSF2",  "CPSF", "CPSF3",  "CPSF",
  "CPSF4",  "CPSF", "WDR33",  "CPSF", "FIP1L1", "CPSF",
  "NUDT21", "CFI", "CPSF6",  "CFI", "CPSF7",  "CFI",
  "PCF11",  "CFII","CLP1",   "CFII",
  "CSTF1",  "CSTF", "CSTF2",  "CSTF", "CSTF3",  "CSTF",
  "CSTF2T", "CSTF",
  "RBBP6",  "Other","SYMPK" ,"Other"                                        
)

# Attach study metadata
cpa_methyl <- Shechter_table %>%
  dplyr::rename(gene_symbol = Gene.Symbol,
         position    = Position.In.Protein,
         reference   = Reference) %>%
  filter(gene_symbol %in% cpa_factors_of_interest) %>%
  # Match reference string to study metadata
  mutate(study_label = NA_character_,
         model       = NA_character_,
         model_type  = NA_character_)

for (i in seq_len(nrow(study_metadata))) {
  matches <- str_detect(cpa_methyl$reference,
                        study_metadata$reference_pattern[i])
  cpa_methyl$study_label[matches] <- study_metadata$study_label[i]
  cpa_methyl$model[matches]       <- study_metadata$model[i]
  cpa_methyl$model_type[matches]  <- study_metadata$model_type[i]
}

# Summarise: n unique sites per gene per study
cpa_summary <- cpa_methyl %>%
  group_by(gene_symbol, study_label, model, model_type) %>%
  summarise(
    n_sites       = n_distinct(position),
    .groups = "drop"
  ) %>%
  left_join(cpa_complex_map, by = "gene_symbol") %>%
  # Ensure all gene x study combinations exist (for empty tiles)
  complete(gene_symbol, study_label,
           fill = list(n_sites = 0))

# ---- Colour palettes --------------------------------------------------------
# Okabe-Ito based — safe for deuteranopia, protanopia, and tritanopia
model_type_colours <- c(
  "Patient tumour tissue" = "#0072B2",   # dark blue
  "Cancer cell line"      = "#D55E00",   # vermillion
  "Immortalised line"          = "#009E73",   # bluish green
  "Database"              = "#999999"    # grey
)

complex_colours <- c(
  "CPSF"  = "#E69F00",   # orange
  "CFI"  = "#56B4E9",   # sky blue
  "CFII" = "#CC79A7",   # rose pink
  "CSTF"  = "#882255",   # dark wine — readable as axis text on white
  "Other" = "#332288"    # dark indigo
)

# Create a numeric model_type rank for sorting
model_type_rank <- c(
  "Patient tumour tissue" = 1,
  "Cancer cell line"      = 2,
  "Immortalised line"          = 3,
  "Database"              = 4
)

# ---- Gene order (fixed, by complex) -----------------------------------------
# New gene order: CFIm, CPSF, CSTF, CFIIm, Other
gene_order <- c(
  "NUDT21", "CPSF6", "CPSF7",                              # CFIm      (1-3)
  "CPSF1", "CPSF2", "CPSF3", "CPSF4", "WDR33", "FIP1L1",   # CPSF (4-9)
  "CSTF1", "CSTF2", "CSTF2T", "CSTF3",                   # CSTF     (10-13)
  "PCF11", "CLP1",                                         # CFIIm    (14-15)
  "RBBP6","SYMPK"                                                  # Other    (16-17)
)

# Vlines after the last gene in each complex
complex_vline_positions <- c(3.5, 9.5, 13.5, 15.5)

# Midpoint of each complex group
complex_label_positions <- c(
  CFI   = 2,
  CPSF   = 6.5,     # midpoint of positions 4-9
  CSTF   = 11.5,  # midpoint of positions 10-13
  CFII  = 14.5,  # midpoint of positions 14-15
  Other  = 16.5
)

# ---- Build combined row label (study + model) --------------------------------
cpa_plot_df <- cpa_summary %>%
  filter(n_sites > 0) %>%
  left_join(
    cpa_methyl %>% distinct(study_label, model, model_type),
    by = "study_label"
  ) %>%
  mutate(
    model      = coalesce(model.y, model.x),
    model_type = coalesce(model_type.y, model_type.x)
  ) %>%
  dplyr::select(-model.x, -model.y, -model_type.x, -model_type.y) %>%
  left_join(cpa_complex_map, by = "gene_symbol") %>%
  mutate(
    row_label       = paste0(study_label, "\n(", model, ")"),
    gene_symbol     = factor(gene_symbol, levels = gene_order),
    model_type      = factor(model_type,
                             levels = c("Patient tumour tissue",
                                        "Cancer cell line",
                                        "Immortalised line",
                                        "Database")),
    model_type_rank = model_type_rank[as.character(model_type)],
    row_label       = fct_reorder(row_label, model_type_rank +
                                    as.numeric(factor(study_label)) * 0.01)
  )

# ---- X-axis label colours by complex ----------------------------------------
x_label_colours <- cpa_complex_map %>%
  filter(gene_symbol %in% gene_order) %>%
  arrange(match(gene_symbol, gene_order)) %>%
  mutate(colour = complex_colours[complex]) %>%
  pull(colour, name = gene_symbol)

# ---- Plot -------------------------------------------------------------------
ggplot(cpa_plot_df, aes(x = gene_symbol, y = row_label)) +
  # Complex separator lines — must come before geom_point
  geom_vline(xintercept = complex_vline_positions,
             linetype = "dashed", colour = "grey70", linewidth = 0.5) +
  geom_point(aes(size = n_sites, colour = model_type), alpha = 0.85) +
  # Complex bracket labels at top using annotate
  annotate("text",
           x     = complex_label_positions,
           y     = Inf,
           label = names(complex_label_positions),
           vjust = 1.5,          # <-- was -0.5, now pulls labels inside the panel
           size = 5, fontface = "bold",
           colour = complex_colours[names(complex_label_positions)])+
  scale_size_continuous(
    name   = "Methylarginine\nsites detected",
    range  = c(2, 10),
    breaks = c(1, 3, 5, 10)
  ) +
  scale_colour_manual(
    name   = "Model system",
    values = model_type_colours
  ) +
  scale_x_discrete(limits = gene_order) +
  coord_cartesian(clip = "off") +   # needed so annotate vjust = -0.5 isn't clipped
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(
      angle  = 45,
      hjust  = 1,
      colour = x_label_colours[gene_order],
      face   = "bold",
      size   = 16
    ),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92"),
    plot.margin = margin(t = 45, r = 10, b = 10, l = 10)
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Arginine methylation of CPA factors across published proteomics studies"
  )

