library(tidyverse)        # purrr, dplyr, readr, ggplot2, etc.

dir_path <- "../../Data/Figure_1_S1_S2_S3/DESeq2_tables/timecourse"

# --------- 2. Find the four result files --------------------------------------
files <- list.files(
  path   = dir_path,
  pattern = "all_genes\\.csv$",        # ends with "all_genes.csv"
  full.names = TRUE
)

if (length(files) == 0) stop("No *all_genes.csv files found in dir_path")


# --------- 3. Read, bind together -------------------------------
pcf11_df <- purrr::map_df(files, function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  # extract the digits right before 'hrs_all_genes.csv'
  tp_hr <- stringr::str_extract(
    basename(f),
    "(?<=_)\\d+(?=hrs_all_genes\\.csv$)"     # e.g. 24 in ..._24hrs_all_genes.csv
  ) |>
    as.numeric()
  
  if (is.na(tp_hr)) stop(paste("Couldn't parse the time‑point in:", f))
  
  dplyr::filter(df, gene_symbol == "PCF11") |>
    dplyr::mutate(time_hr = tp_hr,
                  file    = basename(f))
})

print(pcf11_df)   # optional inspection

# CFIM25
CFIM25_df <- purrr::map_df(files, function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  # extract the digits right before 'hrs_all_genes.csv'
  tp_hr <- stringr::str_extract(
    basename(f),
    "(?<=_)\\d+(?=hrs_all_genes\\.csv$)"     # e.g. 24 in ..._24hrs_all_genes.csv
  ) |>
    as.numeric()
  
  if (is.na(tp_hr)) stop(paste("Couldn't parse the time‑point in:", f))
  
  dplyr::filter(df, gene_symbol == "NUDT21") |>
    dplyr::mutate(time_hr = tp_hr,
                  file    = basename(f))
})

# CSTF2
CSTF2_df <- purrr::map_df(files, function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  # extract the digits right before 'hrs_all_genes.csv'
  tp_hr <- stringr::str_extract(
    basename(f),
    "(?<=_)\\d+(?=hrs_all_genes\\.csv$)"     # e.g. 24 in ..._24hrs_all_genes.csv
  ) |>
    as.numeric()
  
  if (is.na(tp_hr)) stop(paste("Couldn't parse the time‑point in:", f))
  
  dplyr::filter(df, gene_symbol == "CSTF2") |>
    dplyr::mutate(time_hr = tp_hr,
                  file    = basename(f))
})

# CSPF73
CPSF73_df <- purrr::map_df(files, function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  # extract the digits right before 'hrs_all_genes.csv'
  tp_hr <- stringr::str_extract(
    basename(f),
    "(?<=_)\\d+(?=hrs_all_genes\\.csv$)"     # e.g. 24 in ..._24hrs_all_genes.csv
  ) |>
    as.numeric()
  
  if (is.na(tp_hr)) stop(paste("Couldn't parse the time‑point in:", f))
  
  dplyr::filter(df, gene_symbol == "CPSF3") |>
    dplyr::mutate(time_hr = tp_hr,
                  file    = basename(f))
})

# --- combine & order proteins alphabetically for facets ---
all_cpa_df <- bind_rows(
  pcf11_df  %>% mutate(protein = "PCF11"),
  CSTF2_df  %>% mutate(protein = "CSTF2"),
  CPSF73_df %>% mutate(protein = "CPSF73"),   # gene_symbol was "CPSF3"
  CFIM25_df %>% mutate(protein = "CFIM25")    # gene_symbol was "NUDT21"
) %>%
  mutate(
    time_hr  = as.integer(time_hr),
    time_hr_f = factor(time_hr, levels = c(24, 48, 72, 96)),
    protein   = factor(protein, levels = sort(c("PCF11","CSTF2","CPSF73","CFIM25")))
  )

# --- build a colour key: one row per point, one colour per (protein, time) ---
all_cpa_df <- all_cpa_df %>%
  mutate(col_key = interaction(protein, time_hr, sep = "_", drop = TRUE))

# Palettes: PCF11 & CFIM25 = blue shades; CSTF2 = green shades; CPSF73 = orange shades
blue_shades  <- c("24"="#cfe8ff","48"="#9fd0ff","72"="#5ea9ff","96"="#1f78b4")
green_shades <- c("24"="#d7f2d7","48"="#a7e3a7","72"="#5bcf5b","96"="#1b9e1b")
orange_shades<- c("24"="#ffe6cc","48"="#ffcc99","72"="#ff9933","96"="#cc7a00")
purple_shades <-c("24"="#EAD7FF","48"="#C8A8FF","72"="#9965FF","96"="#6A3D9A")  # dark purple

# Update the color map: keep others as before, CFIM25 now uses purple_shades
color_map <- c(
  setNames(blue_shades  [c("24","48","72","96")], paste0("PCF11_",  c(24,48,72,96))),
  setNames(green_shades [c("24","48","72","96")], paste0("CSTF2_",  c(24,48,72,96))),
  setNames(orange_shades[c("24","48","72","96")], paste0("CPSF73_", c(24,48,72,96))),
  setNames(purple_shades[c("24","48","72","96")], paste0("CFIM25_", c(24,48,72,96)))  # ← changed
)

# --- plot: one row of facets; black line; points shaded by time within each protein ---
ggplot(all_cpa_df, aes(x = time_hr, y = log2FoldChange)) +
  geom_line(aes(group = 1), colour = "black") +
  geom_point(aes(colour = col_key), size = 4) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted",
             linewidth = 1, colour = "black") +
  scale_x_continuous(breaks = c(24, 48, 72, 96)) +
  scale_y_continuous(limits = c(-2.75, 2.75)) +
  scale_colour_manual(values = color_map, guide = "none") +
  facet_wrap(~ protein, nrow = 1) +
  labs(
    title = "CPA factor mRNA levels following DMAi",
    x = "DMAi treatment (hours)",
    y = "log2FoldChange"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 20, margin = margin(t = 8)),
    axis.title.y = element_text(size = 20, margin = margin(r = 8)),
    axis.text    = element_text(size = 15),
    strip.text   = element_text(size = 16, face = "bold"),
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
    panel.spacing.x = unit(3, "pt")   # set to unit(0, "pt") if you want borders to touch
  )

