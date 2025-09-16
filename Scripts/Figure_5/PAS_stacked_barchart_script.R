library(tidyverse)
library(stringr)
library(scales) 
library(data.table)

# ── 1. Load .bed files ─────────────────────────────────────────────────────
bed_paths <- list.files(
  path = "../../Data/Figure_4_S7/bed/APA_3seq/200nt_window",
  pattern = "\\.bed$",
  full.names = TRUE,
  ignore.case = TRUE
)

col_titles <- c("chr", "start", "end", "gene_name", "feature_id", "strand")

cfim_mitigation_bed_files <- lapply(
  bed_paths,
  fread,
  sep = "\t",
  header = FALSE,
  col.names = col_titles,
  data.table = FALSE
)

# ── 2. Extract PAS motif from feature_id ───────────────────────────────────
grab_pas <- "(?<=_)\\w+(?=_)"

cfim_mitigation_bed_files <- lapply(
  cfim_mitigation_bed_files,
  function(df) {
    df$PAS <- str_extract(df$feature_id, grab_pas)
    df
  }
)

# ── 3. Define groupings and sample sizes ───────────────────────────────────
idx <- c(3, 4, 5, 6, 7, 8)
condition <- c("DMAi full mitigation", "DMAi full mitigation",
               "DMAi partial mitigation", "DMAi partial mitigation",
               "DMAi no mitigation", "DMAi no mitigation")
location <- c("distal", "proximal",
              "distal", "proximal",
              "distal", "proximal")

row_totals <- tibble(
  idx = idx,
  Condition = condition,
  Location = location,
  total = sapply(cfim_mitigation_bed_files[idx], nrow)
)

# ── 4. Build PAS count + proportion table ──────────────────────────────────
pas_prop_tbl <- map_dfr(
  seq_along(idx),
  function(i) {
    cfim_mitigation_bed_files[[idx[i]]] %>%
      as_tibble() %>%
      dplyr::count(PAS, name = "n") %>%
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
                        c("DMAi no mitigation", "DMAi partial mitigation", "DMAi full mitigation"),
                        " (",
                        row_totals$total[match(
                          c("DMAi no mitigation", "DMAi partial mitigation", "DMAi full mitigation"),
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
    position = "fill", width = 0.7, alpha = 1
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
    strip.text       = element_text(size = 28),
    strip.background = element_rect(fill = "grey90", colour = "black")
  )

#is difference in AAUAAA proportion significant 
# motif_tbl already has one row per PAS_group with counts `n`
aa_counts <- motif_tbl %>%
  filter(Location == "distal pA",
         PAS_group == "AAUAAA",
         Condition %in% c("DMAi full mitigation", "DMAi no mitigation")) %>%
  select(Condition, n) %>%
  pivot_wider(names_from = Condition, values_from = n)   # -> AAUAAA counts

total_counts <- motif_tbl %>%
  filter(Location == "distal pA",
         Condition %in% c("DMAi full mitigation", "DMAi no mitigation")) %>%
  group_by(Condition) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = total)   # -> totals

# pull numbers into vectors
AA_full  <- aa_counts$`DMAi full mitigation`
AA_none  <- aa_counts$`DMAi no mitigation`
tot_full <- total_counts$`DMAi full mitigation`
tot_none <- total_counts$`DMAi no mitigation`

prop.test(c(AA_full, AA_none),
          c(tot_full, tot_none),
          correct = TRUE)            
