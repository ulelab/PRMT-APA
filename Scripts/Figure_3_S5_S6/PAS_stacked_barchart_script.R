library(tidyverse)
library(stringr)
library(scales) 
library(data.table)

# ── 1. Load .bed files ─────────────────────────────────────────────────────
bed_paths <- list.files(
  path       = "../../Data/Figure_3_S5_S6/bed/APA_3seq/200nt_window/DMAi_vs_non_vs_control/",                # use explicit path if not working dir
  pattern    = "*\\.bed$",
  full.names = TRUE,
  ignore.case = TRUE               # drop if you need case–sensitive
)

col_titles <- c("chr", "start", "end", "gene_name", "feature_id", "strand")

bed_files <- lapply(
  bed_paths,
  fread,
  sep = "\t",
  header = FALSE,
  col.names = col_titles,
  data.table = FALSE
)

# ── 2. Extract PAS motif from feature_id ───────────────────────────────────
grab_pas <- "(?<=_)\\w+(?=_)"

bed_files <- lapply(
  bed_files,
  function(df) {
    df$PAS <- str_extract(df$feature_id, grab_pas)
    df
  }
)


# ── 3. Define groupings and sample sizes ───────────────────────────────────
idx <- c(1, 2, 3, 4, 5, 6)
condition <- c("DMAi lengthened", "DMAi lengthened",
               "control", "control",
               "non-DMAi lengthened", "non-DMAi lengthened")
location <- c("proximal", "distal",
              "distal", "proximal",
              "proximal", "distal")

row_totals <- tibble(
  idx = idx,
  Condition = condition,
  Location = location,
  total = sapply(bed_files[idx], nrow)
)

# ── 4. Build PAS count + proportion table ──────────────────────────────────
pas_prop_tbl <- map_dfr(
  seq_along(idx),
  function(i) {
    bed_files[[idx[i]]] %>%
      as_tibble() %>%
      count(PAS, name = "n") %>%
      mutate(
        Condition = condition[i],
        Location = location[i]
      )
  }
) %>%
  group_by(Condition, Location) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# ── 5. Label motifs (separate AAUAAA and AUUAAA) ───────────────────────────
motif_tbl <- pas_prop_tbl %>%
  mutate(
    PAS_group = case_when(
      PAS == "AATAAA" ~ "AAUAAA",
      PAS == "ATTAAA" ~ "AUUAAA",
      PAS == "APA"    ~ "non-canonical",
      PAS == "None"   ~ "None",
      TRUE            ~ PAS
    )
  ) %>%
  group_by(Condition, Location, PAS_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(proportion = n / sum(n))

# ── 6. Join row totals and prepare final labels ────────────────────────────
motif_tbl <- motif_tbl %>%
  left_join(row_totals, by = c("Condition", "Location")) %>%
  mutate(
    Cond_lab = paste0(Condition, " (", total, ")"),
    Cond_lab = factor(Cond_lab,
                      levels = paste0(
                        c("control", "DMAi lengthened", "non-DMAi lengthened"),
                        " (",
                        row_totals$total[match(
                          c("control", "DMAi lengthened", "non-DMAi lengthened"),
                          row_totals$Condition
                        )], ")")),
    Location = factor(Location,
                      levels = c("proximal", "distal"),
                      labels = c("proximal pA", "distal pA"))
  )

# ── 7. Plot stacked bar chart ──────────────────────────────────────────────
ggplot() +
  geom_col(
    data = subset(motif_tbl, Location == "proximal pA"),
    aes(Cond_lab, proportion, fill = PAS_group),
    position = "fill", width = 0.7, alpha = 0.4
  ) +
  geom_col(
    data = subset(motif_tbl, Location == "distal pA"),
    aes(Cond_lab, proportion, fill = PAS_group),
    position = "fill", width = 0.7, alpha = 1
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  facet_wrap(~Location, nrow = 1) +
  scale_fill_manual(
    values = c(
      "AAUAAA"        = "#377eb8",
      "AUUAAA"        = "#4daf4a",
      "non-canonical" = "#e41a1c",
      "None"          = "#999999"
    )
  ) +
  labs(x = NULL,
       y = "Percentage of PAS motifs",
       fill = "PAS class") +
  theme_classic(base_size = 14) +
  theme(
    text             = element_text(size = 18),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(size = 18),
    strip.background = element_rect(fill = "grey90", colour = "black")
  )
