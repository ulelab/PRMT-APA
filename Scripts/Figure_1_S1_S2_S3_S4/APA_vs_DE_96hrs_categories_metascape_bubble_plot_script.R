#load libraries
library(tidyverse)

# 1) Read in file and retain relevant columns
metascape_data <- read.csv('Data/Figure_1_S1_S2_S3_S4/APA_vs_DE_96hrs_filtered_collated_metascape_results.csv.csv')

metascape_data = metascape_data[,1:6]

# 2) Unite category column
metascape_data <- metascape_data %>%
  unite("GO Category", GO, Description, sep = ":")

# 3) Create -log10(p) column
metascape_data <- metascape_data %>%
  mutate(`-log10(p)` = -1 * Log10.P.) %>%
  dplyr::select(-Log10.P.)

#reorder factor
metascape_data <- metascape_data %>%
  mutate(
    Group = recode(
      Group,
      "down_and_lengthened"      = "3' UTR lengthening + gene downregulated",
      "up_and_lengthened"        = "3' UTR lengthening + gene upregulated",
      "lengthened_no_change"     = "3' UTR lengthening + no expression change",
      "down_not_lengthened"      = "no 3' UTR lengthening + gene downregulated",
      "up_not_lengthened"        = "no 3' UTR lengthening + gene upregulated"
    ),
    # Make `Group` a factor in the order you want them to appear
    Group = factor(
      Group,
      levels = c(
        "3' UTR lengthening + no expression change",
        "3' UTR lengthening + gene upregulated",
        "no 3' UTR lengthening + gene upregulated",
        "3' UTR lengthening + gene downregulated",
        "no 3' UTR lengthening + gene downregulated"
      )
    )
  )

# 4) Calculate the maximum -log10(p) per GO Category, create GroupCategory variable
metascape_data <- metascape_data %>%
  group_by(Group) %>%
  # rank: higher -log10(p) => higher rank number
  mutate(rank_within_group = rank(`-log10(p)`, ties.method = "first")) %>%
  ungroup() %>%
  # paste Group name + GO Category into one label
  mutate(GroupCategory = paste(Group, `GO Category`, sep = " | ")) %>%
  # arrange by Group (in factor order) and then rank
  arrange(Group, rank_within_group) %>%
  # turn GroupCategory into a factor that follows this arranged order
  mutate(GroupCategory = factor(GroupCategory, levels = unique(GroupCategory)))

# 5) Define a custom function to split long labels
split_label <- function(x, width = 100) {
  # strwrap() returns a list of vectors, each vector is 
  # the wrapped lines for one input string
  lines_list <- strwrap(x, width = width, simplify = FALSE)
  # For each element in lines_list, collapse lines with "\n"
  sapply(lines_list, function(lines) paste(lines, collapse = "\n"))
}

# 6) Define a named vector of colors mapped to your Group column
my_colors <- c(
  "3' UTR lengthening + gene upregulated"       = "#670303",
  "3' UTR lengthening + gene downregulated"     = "#f98a87",
  "3' UTR lengthening + no expression change"   = "#009E73",
  "no 3' UTR lengthening + gene upregulated"    = "#004278",
  "no 3' UTR lengthening + gene downregulated"  = "#5dbfed"
)

#create boundaries for divider lines
df_boundaries <- metascape_data %>%
  distinct(Group, GroupCategory) %>%
  group_by(Group) %>%
  summarize(boundary = max(as.numeric(GroupCategory))) %>%
  ungroup() %>%
  arrange(boundary)

# If you don't want a line after the last group, remove the final row:
df_boundaries <- df_boundaries[-nrow(df_boundaries), ]

# 7) Plot using color = Group
ggplot(metascape_data, aes(
  x = `-log10(p)`,
  y = GroupCategory,
  color = Group,
  size  = Count
)) +
  geom_point() +
  scale_color_manual(values = my_colors) +
  scale_size(range = c(4, 11)) +
  labs(x = "-log10(p)", y = "GO Category") +
  geom_hline(
    data = df_boundaries,
    aes(yintercept = boundary + 0.5),
    color = "grey40", linetype = "dotted"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title  = element_text(size = 16),
    axis.text.y = element_text(size = 11),
    legend.title= element_text(size = 13),
    legend.text = element_text(size = 13)
  ) +
  scale_y_discrete(labels = function(x) {
    # 1) Remove everything up to and including ' | '
    short_label <- sub(".*\\| ", "", x)
    # 2) Wrap to 50 characters
    split_label(short_label, 80)
  })

