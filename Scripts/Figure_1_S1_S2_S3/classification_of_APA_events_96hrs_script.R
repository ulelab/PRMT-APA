# Load required libraries
library(tidyverse)   # For data manipulation and visualization (includes dplyr, ggplot2)
library(data.table)  # For efficient data manipulation
library(RColorBrewer)# For color palettes in plots
library(ggh4x)       # Extensions for ggplot2
library(stringr)     # For string manipulation

# ----------------------------
# Set Working Directories and File Paths
# ----------------------------

# Update these paths to match your data locations
polyA_csv_dir <- "/path/to/your/CSV_files"
anno_bed_dir <- "/path/to/your/annotated_BED_files"
UTR_info_path <- "/path/to/your/UTR_regions.bed"
splice_site_info_path <- "/path/to/your/splice_regions.bed"
output_dir <- "/path/to/your/output_directory"

# For demonstration purposes, replace the above paths with your actual paths
# For example:
polyA_csv_dir <- "/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/dedup/github/with_intronic/CSV_files"
anno_bed_dir <- "/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/dedup/github/with_intronic/TUTR_vs_ALE"
UTR_info_path <- "/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/no_dedup/two_step_test/TUTR_vs_ALE/UTR_regions_gencode_v45.bed"
splice_site_info_path <- "/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/no_dedup/two_step_test/TUTR_vs_ALE/splice_regions_gencode_v45.bed"
output_dir <- "/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/dedup/github/with_intronic/TUTR_vs_ALE"

# Ensure the output directory exists
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the experimental conditions
conditions <- c("ADMAi", "SDMAi_ADMAi","SDMAi")  # Update this vector as per your conditions

# ----------------------------
# Read and Prepare PolyA Site Data
# ----------------------------

# Set working directory to where your PolyA site CSV files are located
setwd(polyA_csv_dir)

# List all CSV files containing 'all' in the current working directory
polyA_files <- list.files(pattern = "all.*\\.csv$", full.names = TRUE)

# Read each file into a list of data frames
polyA_data_list <- lapply(polyA_files, read.csv)

# Tidy up 'ensembl_ID' column by removing version numbers or any suffix after the dot
polyA_data_list <- lapply(polyA_data_list, function(df){
  df$ensembl_ID <- sub("\\..*$", "", df$ensembl_ID)
  return(df)
})

# Calculate mean usage of each PolyA site across control and treatment conditions
polyA_data_list <- lapply(polyA_data_list, function(df){
  df$mean_usage <- (df$pA_usage_control + df$pA_usage_treatment) / 2 
  return(df)
})

# Function to process PolyA site data and select top 2 sites per gene
process_gene_groups <- function(df) {
  # Add a flag for significant PolyA sites based on adjusted p-value and change in usage
  df <- df %>%
    mutate(is_significant = ifelse(!is.na(twostep_feature_padj) & 
                                     twostep_feature_padj < 0.05 & 
                                     abs(change_in_usage) > 0.1, TRUE, FALSE))
  
  # Count the number of significant PolyA sites per gene
  df <- df %>%
    group_by(gene_name) %>%
    mutate(sig_PAS_count = sum(is_significant, na.rm = TRUE)) %>%
    ungroup()
  
  # For genes with at least two significant PAS
  df_sig <- df %>%
    dplyr::filter(sig_PAS_count >= 2 & is_significant) %>%  # Keep only significant PAS
    group_by(gene_name) %>%
    arrange(desc(abs(change_in_usage))) %>%  # Order by absolute change in usage
    slice_head(n = 2) %>%  # Select top 2 PAS per gene
    ungroup()
  
  # For genes with fewer than two significant PAS
  df_non_sig <- df %>%
    dplyr::filter(sig_PAS_count < 2) %>%
    group_by(gene_name) %>%
    arrange(desc(mean_usage)) %>%  # Order by mean usage
    slice_head(n = 2) %>%  # Select top 2 PAS per gene
    ungroup()
  
  # Combine the two data frames
  df_result <- bind_rows(df_sig, df_non_sig) %>% dplyr::select(-is_significant)
  
  return(df_result)
}

# Apply the function to each data frame in the list
polyA_top2_list <- lapply(polyA_data_list, process_gene_groups)

# ----------------------------
# Read and Prepare Annotated BED Files
# ----------------------------

# Set working directory to where your annotated BED files are located
setwd(anno_bed_dir)

# Get a list of all BED files in the directory
bed_files <- list.files(pattern = '\\.bed$', full.names = TRUE)

# Read each BED file into a list of data frames
anno_data_list <- lapply(bed_files, function(file) {
  read.delim(file, header = FALSE, stringsAsFactors = FALSE,
             col.names = c("chr", "start", "end", "col4", "gene_name", "strand",
                           "gtf_chr", "database", "region", "gtf_start",
                           "gtf_end", "col12", "gtf_strand", "col14", "gtf_annotation"))
})

# Separate 'gtf_annotation' column to extract 'ensembl_ID' and 'extra_annotation'
anno_data_list <- lapply(anno_data_list, function(df) {
  df %>%
    separate(gtf_annotation, into = c("ensembl_ID", "extra_annotation"), sep = ";", extra = "merge", fill = "right") %>%
    mutate(
      ensembl_ID = gsub("^gene_id ", "", ensembl_ID),  # Remove 'gene_id ' prefix
      ensembl_ID = sub("\\..*", "", ensembl_ID)        # Remove version numbers
    )
})

# Ensure that 'polyA_top2_list' and 'anno_data_list' are in the same order
# Join the processed PolyA data with the annotated sites
anno_data_list <- Map(function(df1, df2) {
  df1 %>%                                   # polyA_top2 element
    dplyr::select(-ensembl_ID) %>%          # drop its ensembl_id first
    dplyr::left_join(                       # join to anno element
      df2, 
      by = c("chr", "start", "end", "strand", "gene_name")
    )
},
polyA_top2_list,
anno_data_list
)

# Further separate 'extra_annotation' to extract 'transcript_annotation' and 'further_annotation'
anno_data_list <- lapply(anno_data_list, function(df) {
  df %>%
    mutate(
      transcript_annotation = if_else(
        grepl("transcript_id", extra_annotation),
        sub(";.*", "", extra_annotation),
        extra_annotation
      ),
      further_annotation = if_else(
        grepl("transcript_id", extra_annotation),
        sub(".*?;", "", extra_annotation),
        NA_character_
      )
    ) %>%
    dplyr::select(-extra_annotation)
})

# If UTR annotation is present in each 'transcript_annotation' group, only keep rows with 'UTR' in 'region'
anno_data_UTR_filtered <- lapply(anno_data_list, function(df) {
  df %>%
    group_by(chr, start, end, transcript_annotation) %>%
    dplyr::filter(!any(grepl("UTR", region)) | grepl("UTR", region)) %>%
    ungroup()
})

# Remove 'gene' annotations if there are other annotations for that PolyA site
anno_data_UTR_filtered <- lapply(anno_data_UTR_filtered, function(df) {
  df %>%
    group_by(chr, start, end) %>%
    mutate(has_extra_anno = any(region != 'gene')) %>%
    dplyr::filter(!(region == 'gene' & has_extra_anno)) %>%
    dplyr::select(-has_extra_anno) %>%
    ungroup()
})

# Remove 'transcript' annotations if there are other annotations for that PolyA site
anno_data_UTR_filtered <- lapply(anno_data_UTR_filtered, function(df) {
  df %>%
    group_by(chr, start, end, ensembl_ID) %>%
    mutate(has_extra_anno = any(region != 'transcript')) %>%
    dplyr::filter(!(region == 'transcript' & has_extra_anno)) %>%
    dplyr::select(-has_extra_anno) %>%
    ungroup()
})

# Read in UTR location file
UTR_info <- read.delim(UTR_info_path, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("UTR_chr", "UTR_start", "UTR_end", "UTR_strand", "ensembl_ID", "UTR_class", "transcript_annotation"))

# Keep only 3' UTR information
UTR3_info <- UTR_info %>% 
  dplyr::filter(UTR_class == '3UTR')

# Prepare UTR3_info for merging
UTR3_info$UTR_chr <- paste0("chr", UTR3_info$UTR_chr)
UTR3_info$ensembl_ID <- sub("\\..*$", "", UTR3_info$ensembl_ID)
UTR3_info$transcript_annotation <- paste(" transcript_id ", UTR3_info$transcript_annotation, sep = "")

# Merge annotated sites with UTR start and end information
anno_data_with_UTR <- lapply(anno_data_UTR_filtered, function(df) {
  dplyr::left_join(df, UTR3_info, by = c("ensembl_ID", "transcript_annotation"))
})

# Split data frames into positive and negative strands
anno_data_pos <- lapply(anno_data_with_UTR, function(df) {
  df %>% dplyr::filter(strand == '+')
})
anno_data_neg <- lapply(anno_data_with_UTR, function(df) {
  df %>% dplyr::filter(strand == '-')
})

#read in gtf annotation file
gencode_v45 = read.delim('/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/own_atlas/TUTR_vs_ALE/gencode.v45.annotation.gtf', stringsAsFactors = FALSE,
                         col.names = c("chr", "database", "region", "gtf_start", "gtf_end", "col6",
                                       "strand", "col8", "gtf_annotation"))

#get annotation column in same structure as other files
gencode_v45 = gencode_v45 %>%
  separate(gtf_annotation, into = c("ensembl_ID", "extra_annotation"), sep = ";", extra = "merge", fill = "right") %>%
  mutate(
    ensembl_ID = gsub("^gene_id ", "", ensembl_ID),  # Remove 'gene_id ' prefix
    ensembl_ID = sub("\\..*", "", ensembl_ID)        # Remove version numbers
  )

gencode_v45 = gencode_v45 %>%
  mutate(
    transcript_annotation = if_else(
      grepl("transcript_id", extra_annotation),
      sub(";.*", "", extra_annotation),
      extra_annotation
    ),
    further_annotation = if_else(
      grepl("transcript_id", extra_annotation),
      sub(".*?;", "", extra_annotation),
      NA_character_
    )
  ) %>%
  dplyr::select(-extra_annotation)

#find 3' end of all transcripts (+ve strand)
gencode_max_gtf_pos <- gencode_v45 %>%
  group_by(transcript_annotation) %>%
  summarise(max_gtf_end = max(gtf_end, na.rm = TRUE)) %>%
  ungroup()

#find 3' end of all transcripts (-ve strand)
gencode_min_gtf_neg <- gencode_v45 %>%
  group_by(transcript_annotation) %>%
  summarise(min_gtf_start = min(gtf_start, na.rm = TRUE)) %>%
  ungroup()

# For positive strand, keep annotations where transcript end is near PolyA site and that are high TSL
boundary_plus_TSL_filter_pos <- function(df, gencode_max_gtf) {
  
  # 1) Join the input df to the precomputed 'max_gtf_end' from gencode_v45
  df_joined <- df %>%
    left_join(gencode_max_gtf, by = "transcript_annotation")
  
  # 2) Define boundary_value:
  #    If UTR_end != NA, use UTR_end
  #    Otherwise, use max_gtf_end
  df_joined <- df_joined %>%
    mutate(
      boundary_value = if_else(is.na(UTR_end), max_gtf_end, UTR_end),
      keep_row = boundary_value > (end - 100) & boundary_value < (end + 30)
    ) %>%
    filter(keep_row) %>%
    select(-boundary_value, -keep_row, -max_gtf_end)
  
  # 3) TSL filtering on the rows that passed the boundary check
  df_in_boundary <- df_joined %>%
    mutate(further_annotation = coalesce(further_annotation, "None")) %>%
    group_by(chr, start, end, gene_name) %>%
    filter(
      if (any(str_detect(further_annotation, "transcript_support_level 1"))) {
        str_detect(further_annotation, "transcript_support_level 1") | 
          str_detect(further_annotation, "transcript_support_level 2")
      } else if (any(str_detect(further_annotation, "transcript_support_level 2"))) {
        str_detect(further_annotation, "transcript_support_level 2")
      } else if (any(str_detect(further_annotation, "transcript_support_level 3"))) {
        str_detect(further_annotation, "transcript_support_level 3")
      } else {
        TRUE
      }
    ) %>%
    ungroup()
  
  # 4) Mark these rows as near the UTR boundary
  df_in_boundary <- df_in_boundary %>%
    mutate(near_UTR_boundary = "YES")
  
  df_in_boundary
}

# run filter
anno_data_pos_filtered <- lapply(anno_data_pos, function(df_elem) {
  boundary_plus_TSL_filter_pos(df_elem, gencode_max_gtf_pos)
})

# For each condition i, build a table of all transcript_annotations that passed for each gene
passed_annos_pos <- lapply(anno_data_pos_filtered, function(df) {
  df %>%
    group_by(gene_name) %>%
    summarise(passed_annos = list(unique(transcript_annotation))) %>%
    ungroup()
})

#find the sites that failed the boundary test
lost_sites_pos <- Map(function(original_df, filtered_df){
  lost_df = anti_join(original_df, filtered_df, by = c("chr","start","end"))
  lost_df = lost_df %>%
    mutate(near_UTR_boundary = "NO")
  lost_df
}, anno_data_pos, anno_data_pos_filtered)

#find the starting point of the final exon all transcripts (+ve strand)
gencode_max_gtf_start_pos <- gencode_v45 %>%
  group_by(transcript_annotation) %>%
  summarise(max_gtf_start = max(gtf_start, na.rm = TRUE)) %>%
  ungroup()

#find the starting point of the final exon all transcripts (-ve strand)
gencode_min_gtf_end_neg <- gencode_v45 %>%
  group_by(transcript_annotation) %>%
  summarise(min_gtf_end = min(gtf_end, na.rm = TRUE)) %>%
  ungroup()


#rescue these sites by matching the annotations to the other polyA site of that gene, or fallback use TSL and distance
rescue_lost_sites_pos <- function(lost_df,
                                  passed_anno_df,
                                  gencode_max_gtf_start_pos) {
  ## 1. Flatten the passed‑annotation table
  passed_flat <- passed_anno_df %>%
    tidyr::unnest(cols = passed_annos, names_repair = "unique") %>%
    dplyr::rename(passed_annotation = tidyselect::matches("^passed_annos"))
  
  ## 2. Split lost sites by polyA‑site group
  group_list <- lost_df %>%
    group_by(chr, start, end, gene_name) %>%
    group_split()
  
  ## 3. Process each group
  rescue_list <- lapply(group_list, function(df_group) {
    
    ## --- Attempt adoption ---
    adopt_df <- df_group %>%
      left_join(passed_flat, by = "gene_name") %>%
      filter(transcript_annotation == passed_annotation) %>%
      select(-passed_annotation)
    
    if (nrow(adopt_df) > 0) {
      return(adopt_df)          # → adoption succeeded
    }
    
    ## --- Fallback: TSL filter then nearest‑end selection ---
    fallback_df <- df_group %>%
      mutate(further_annotation = coalesce(further_annotation, "None")) %>%
      group_by(chr, start, end, gene_name) %>%
      filter(
        if (any(str_detect(further_annotation, "transcript_support_level 1"))) {
          str_detect(further_annotation, "transcript_support_level 1") |
            str_detect(further_annotation, "transcript_support_level 2")
        } else if (any(str_detect(further_annotation, "transcript_support_level 2"))) {
          str_detect(further_annotation, "transcript_support_level 2")
        } else if (any(str_detect(further_annotation, "transcript_support_level 3"))) {
          str_detect(further_annotation, "transcript_support_level 3")
        } else {
          TRUE
        }
      ) %>%
      ungroup() %>%
      mutate(dist_to_end = if_else(is.na(UTR_end),
                                   abs(gtf_end - end),
                                   abs(UTR_end - end))) %>%
      group_by(chr, start, end, gene_name) %>%
      filter(dist_to_end == min(dist_to_end)) %>%
      ungroup() %>%
      select(-dist_to_end)
    
    ## >>> NEW: update near_UTR_boundary using gencode_max_gtf_start_pos <<<
    fallback_df <- fallback_df %>%
      left_join(gencode_max_gtf_start_pos,
                by = "transcript_annotation") %>%
      mutate(
        near_UTR_boundary = if_else(
          !is.na(max_gtf_start) & start > max_gtf_start,
          "YES",
          near_UTR_boundary        # leave as‑is (usually "NO")
        )
      ) %>%
      select(-max_gtf_start)      # drop helper column
    
    fallback_df                  # return the rescued rows
  })
  
  ## 4. Bind all rescued groups together
  bind_rows(rescue_list)
}

#apply filter
lost_sites_pos_filtered <- Map(
  rescue_lost_sites_pos,
  lost_sites_pos,          # list of lost‑site dfs
  passed_annos_pos,        # list of passed‑annotation dfs
  MoreArgs = list(gencode_max_gtf_start_pos = gencode_max_gtf_start_pos)
)

#combine lost with main table
anno_data_pos_combined <- Map(rbind, anno_data_pos_filtered, lost_sites_pos_filtered)

# if we have multiple gene annotations per polyA site, keep the one closest to the end of the gtf annotation (as long as these are UTRs)
anno_data_pos_combined = lapply(anno_data_pos_combined, function(df) {
  df %>% 
    group_by(feature_id) %>%                                   # ── 1. one pA site
    group_modify(~ {
      ## Keep everything if at least one line isn’t a UTR
      if (!all(.x$region == "UTR")) return(.x)
      
      ## 2.  work at transcript (ensembl_ID) level
      winners <- .x %>% 
        group_by(ensembl_ID) %>%                               # one distance per tx
        summarise(dist = min(abs(end - gtf_end)),
                  .groups = "drop") %>% 
        filter(dist == min(dist)) %>%                          # 3. keep ties
        pull(ensembl_ID)
      
      ## 4. retain *all* rows belonging to the winning transcript(s)
      .x %>% filter(ensembl_ID %in% winners)
    }) %>% 
    ungroup()
})

##same for negative sites
boundary_plus_TSL_filter_neg <- function(df, gencode_min_gtf) {
  
  # 1) Join the input df to the precomputed 'min_gtf_start' from gencode_v45
  df_joined <- df %>%
    left_join(gencode_min_gtf, by = "transcript_annotation")
  
  # 2) Define boundary_value:
  #    If UTR_start != NA, use UTR_start
  #    Otherwise, use min_gtf_start
  df_joined <- df_joined %>%
    mutate(
      boundary_value = if_else(is.na(UTR_start), min_gtf_start, UTR_start),
      keep_row = boundary_value < (start + 100) & boundary_value > (start - 30)
    ) %>%
    filter(keep_row) %>%
    select(-boundary_value, -keep_row, -min_gtf_start)
  
  # 3) TSL filtering on the rows that passed the boundary check
  df_in_boundary <- df_joined %>%
    mutate(further_annotation = coalesce(further_annotation, "None")) %>%
    group_by(chr, start, end, gene_name) %>%
    filter(
      if (any(str_detect(further_annotation, "transcript_support_level 1"))) {
        str_detect(further_annotation, "transcript_support_level 1") | 
          str_detect(further_annotation, "transcript_support_level 2")
      } else if (any(str_detect(further_annotation, "transcript_support_level 2"))) {
        str_detect(further_annotation, "transcript_support_level 2")
      } else if (any(str_detect(further_annotation, "transcript_support_level 3"))) {
        str_detect(further_annotation, "transcript_support_level 3")
      } else {
        TRUE
      }
    ) %>%
    ungroup()
  
  # 4) Mark these rows as near the UTR boundary
  df_in_boundary <- df_in_boundary %>%
    mutate(near_UTR_boundary = "YES")
  
  df_in_boundary
}

# run same functions for negative strand
anno_data_neg_filtered <- lapply(anno_data_neg, function(df_elem) {
  boundary_plus_TSL_filter_neg(df_elem, gencode_min_gtf_neg)
})

passed_annos_neg <- lapply(anno_data_neg_filtered, function(df) {
  df %>%
    group_by(gene_name) %>%
    summarise(passed_annos = list(unique(transcript_annotation))) %>%
    ungroup()
})

lost_sites_neg <- Map(function(original_df, filtered_df){
  lost_df = anti_join(original_df, filtered_df, by = c("chr","start","end"))
  lost_df = lost_df %>%
    mutate(near_UTR_boundary = "NO")
  lost_df
}, anno_data_neg, anno_data_neg_filtered)


rescue_lost_sites_neg <- function(lost_df,
                                  passed_anno_df,
                                  gencode_min_gtf_end_neg) {
  ## ‑‑‑ 1.  flatten the passed‑annotation list ‑‑‑
  passed_flat <- passed_anno_df %>%
    tidyr::unnest(cols = passed_annos, names_repair = "unique") %>%
    dplyr::rename(passed_annotation = tidyselect::matches("^passed_annos"))
  
  ## ‑‑‑ 2.  split lost sites by polyA‑site group ‑‑‑
  group_list <- lost_df %>%
    group_by(chr, start, end, gene_name) %>%
    group_split()
  
  ## ‑‑‑ 3.  process each group ‑‑‑
  rescue_list <- lapply(group_list, function(df_group) {
    
    ## ──‑ Attempt adoption
    adopt_df <- df_group %>%
      left_join(passed_flat, by = "gene_name") %>%
      filter(transcript_annotation == passed_annotation) %>%
      select(-passed_annotation)
    
    if (nrow(adopt_df) > 0) {
      return(adopt_df)              # adoption succeeded
    }
    
    ## ──‑ Fallback: TSL filter → nearest‑end row
    fallback_df <- df_group %>%
      mutate(further_annotation = coalesce(further_annotation, "None")) %>%
      group_by(chr, start, end, gene_name) %>%
      filter(
        if (any(str_detect(further_annotation, "transcript_support_level 1"))) {
          str_detect(further_annotation, "transcript_support_level 1") |
            str_detect(further_annotation, "transcript_support_level 2")
        } else if (any(str_detect(further_annotation, "transcript_support_level 2"))) {
          str_detect(further_annotation, "transcript_support_level 2")
        } else if (any(str_detect(further_annotation, "transcript_support_level 3"))) {
          str_detect(further_annotation, "transcript_support_level 3")
        } else {
          TRUE
        }
      ) %>%
      ungroup() %>%
      ## nearest 3′ end for negative strand
      mutate(dist_to_end = if_else(is.na(UTR_start),
                                   abs(gtf_start - start),
                                   abs(UTR_start - start))) %>%
      group_by(chr, start, end, gene_name) %>%
      filter(dist_to_end == min(dist_to_end)) %>%
      ungroup() %>%
      select(-dist_to_end)
    
    ## ──‑ NEW: promote to boundary if end < min_gtf_end
    fallback_df <- fallback_df %>%
      left_join(gencode_min_gtf_end_neg, by = "transcript_annotation") %>%
      mutate(
        near_UTR_boundary = if_else(
          !is.na(min_gtf_end) & end < min_gtf_end,
          "YES",
          near_UTR_boundary          # leave as is (usually "NO")
        )
      ) %>%
      select(-min_gtf_end)
    
    fallback_df
  })
  
  ## ‑‑‑ 4.  bind all rescued groups together ‑‑‑
  bind_rows(rescue_list)
}

lost_sites_neg_filtered <- Map(
  rescue_lost_sites_neg,
  lost_sites_neg,           # list of lost‑site data frames (−ve strand)
  passed_annos_neg,         # list of passed‑annotation data frames
  MoreArgs = list(gencode_min_gtf_end_neg = gencode_min_gtf_end_neg)
)

anno_data_neg_combined <- Map(rbind, anno_data_neg_filtered, lost_sites_neg_filtered)

# if we have multiple gene annotations per polyA site, keep the one closest to the end of the gtf annotation (as long as these are UTRs)
anno_data_neg_combined = lapply(anno_data_neg_combined, function(df) {
  df %>% 
    group_by(feature_id) %>%                                   # ── 1. one pA site
    group_modify(~ {
      ## Keep everything if at least one line isn’t a UTR
      if (!all(.x$region == "UTR")) return(.x)
      
      ## 2.  work at transcript (ensembl_ID) level
      winners <- .x %>% 
        group_by(ensembl_ID) %>%                               # one distance per tx
        summarise(dist = min(abs(start - gtf_start)),
                  .groups = "drop") %>% 
        filter(dist == min(dist)) %>%                          # 3. keep ties
        pull(ensembl_ID)
      
      ## 4. retain *all* rows belonging to the winning transcript(s)
      .x %>% filter(ensembl_ID %in% winners)
    }) %>% 
    ungroup()
})

# combine positive and negative strand data
anno_data_combined <- Map(rbind, anno_data_pos_combined, anno_data_neg_combined)

# Read in splice site location file
splice_site_info <- read.delim(splice_site_info_path, header = FALSE, stringsAsFactors = FALSE,
                               col.names = c("splice_chr", "splice_start", "splice_end", "coordinates", "splice_score", "splice_strand", "splice_site_type", "transcript_annotation", "ensembl_ID"))

# Keep only final splice site for each transcript
final_splice_site_info <- splice_site_info %>%
  group_by(transcript_annotation) %>%
  dplyr::filter(ifelse(splice_strand == '+', splice_start == max(splice_start), splice_end == min(splice_end))) %>%
  ungroup()

# Prepare splice site info for merging
final_splice_site_info$splice_chr <- paste0("chr", final_splice_site_info$splice_chr)
final_splice_site_info$ensembl_ID <- sub("\\..*$", "", final_splice_site_info$ensembl_ID)
final_splice_site_info$transcript_annotation <- paste(" transcript_id ", final_splice_site_info$transcript_annotation, sep = "")

# Merge annotated sites with splice site information
sites_with_splice_info_list <- lapply(anno_data_combined, function(df) {
  dplyr::left_join(df, final_splice_site_info, by = c("ensembl_ID", "transcript_annotation"))
})

#because of new annotation, remove single pA genes, and if any are three, select two most changing/most used
## helper: operate on ONE ensembl_ID ---------------------------------------
pick_two_features <- function(gene_df) {
  # 1. summarise per PAS
  feature_tbl <- gene_df %>% 
    group_by(feature_id) %>% 
    summarise(
      sig        = any(!is.na(twostep_feature_padj) &
                         twostep_feature_padj < 0.05 &
                         abs(change_in_usage) > 0.1),
      max_change = max(abs(change_in_usage), na.rm = TRUE),
      mean_usage = mean(mean_usage,              na.rm = TRUE),
      .groups    = "drop"
    )
  
  # 2. choose which two PAS to keep
  keep_feats <- if (sum(feature_tbl$sig) >= 2) {
    feature_tbl %>% 
      filter(sig) %>% 
      slice_max(max_change, n = 2, with_ties = FALSE) %>% 
      pull(feature_id)
  } else {
    feature_tbl %>% 
      slice_max(mean_usage, n = 2, with_ties = FALSE) %>% 
      pull(feature_id)
  }
  
  # 3. return ALL original rows that belong to those PAS
  gene_df %>% filter(feature_id %in% keep_feats)
}

## main wrapper: apply to every data frame in your list --------------------
sites_with_splice_info_list <- lapply(
  sites_with_splice_info_list,
  function(df) {
    df %>% 
      ## ── ❶ discard genes whose Δusage is entirely NA ──────────────────
      group_by(ensembl_ID) %>% 
      filter(any(!is.na(change_in_usage))) %>% 
      ## ── ❷ discard genes with only one PAS ───────────────────────────
      filter(n_distinct(feature_id) > 1) %>% 
      ungroup() %>% 
      ## ── ❸ split → process → bind back ───────────────────────────────
      group_split(ensembl_ID, .keep = TRUE) %>% 
      map_dfr(pick_two_features)
  }
)


# Add 'UTR_type' column based on criteria
sites_with_UTR_type_list <- lapply(sites_with_splice_info_list, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    mutate(
      UTR_type = case_when(
        n_distinct(transcript_annotation) == 1 ~ 'TUTR',
        all(is.na(splice_start)) ~ 'TUTR',
        TRUE ~ 'UNKNOWN'
      )
    ) %>%
    ungroup()
})

# Separate 'TUTR' sites
TUTR_sites_list <- lapply(sites_with_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'TUTR')
})

# Process 'UNKNOWN' sites
unknown_sites_list <- lapply(sites_with_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'UNKNOWN')
})

# If there is only one terminal splice site per gene, classify as 'TUTR'
sites_with_UTR_type_list <- lapply(unknown_sites_list, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    mutate(
      UTR_type = case_when(
        n_distinct(splice_start) == 1 ~ 'TUTR',
        TRUE ~ 'UNKNOWN'
      )
    ) %>%
    ungroup()
})

# Add newly identified 'TUTR' sites
more_TUTR_sites_list <- lapply(sites_with_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'TUTR')
})

# Combine 'TUTR' sites
TUTR_sites_list <- Map(rbind, TUTR_sites_list, more_TUTR_sites_list)

# Identify 'iAPA' sites within 'TUTR' sites
TUTR_sites_list <- lapply(TUTR_sites_list, function(df) {
  df %>%
    group_by(ensembl_ID) %>%
    mutate(
      UTR_type = case_when(
        strand == "+" & splice_start >= min(start, na.rm = TRUE) & splice_start <= max(start, na.rm = TRUE) ~ "iAPA",
        strand == "-" & splice_end >= min(end, na.rm = TRUE) & splice_end <= max(end, na.rm = TRUE) ~ "iAPA",
        TRUE ~ UTR_type
      )
    ) %>%
    ungroup()
})

# Separate 'iAPA' sites
iAPA_sites_list <- lapply(TUTR_sites_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'iAPA')
})

# Update 'TUTR' sites
TUTR_sites_list <- lapply(TUTR_sites_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'TUTR')
})

# Process remaining 'UNKNOWN' sites
remaining_unknown_sites_list <- lapply(sites_with_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'UNKNOWN')
})

#find more iAPA sites based on those pA sites that aren't located at the 3' end boundary of their assigned annotation and before the final splice site
remaining_unknown_sites_list = lapply(remaining_unknown_sites_list,function(df){
  df <- df %>%
  # 1) Update UTR_type to iAPA where near_UTR_boundary == "NO" AND
  #    ((strand == "+" AND end < splice_end) OR (strand == "-" AND start > splice_start))
  mutate(UTR_type = case_when(
    near_UTR_boundary == "NO" & strand == "+" & end < splice_end ~ "iAPA",
    near_UTR_boundary == "NO" & strand == "-" & start > splice_start ~ "iAPA",
    TRUE ~ UTR_type
  )) %>%
  # 2) If ANY row in a gene_name group is iAPA, make ALL rows in that gene_name group iAPA.
  group_by(ensembl_ID) %>%
  mutate(
    # Use a plain 'if' to set the entire group’s column if *any* row is iAPA
    UTR_type = if ( any(UTR_type == "iAPA") ) {
      "iAPA"
    } else {
      UTR_type
    }
  ) %>%
  ungroup() 
})

# Separate 'iAPA' sites
more_iAPA_sites_list <- lapply(remaining_unknown_sites_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'iAPA')
})

# Combine 'iAPA' sites
iAPA_sites_list <- Map(rbind, iAPA_sites_list, more_iAPA_sites_list)

#if there are potential TUTR annotations in these groups, we should switch to MIXED. first ID the potential TUTR
iAPA_sites_list = lapply(iAPA_sites_list,function(df) {
  df %>% 
    group_by(transcript_annotation) %>% 
    mutate(
      UTR_type = case_when(
        # ───────── + strand rule ────────────────────────────────────────────
        dplyr::first(strand) == "+" &
          n_distinct(end  [end   > dplyr::first(splice_end)]) >= 2 ~ "TUTR",
        
        # ───────── – strand rule ────────────────────────────────────────────
        dplyr::first(strand) == "-" &
          n_distinct(start[start < dplyr::first(splice_start)]) >= 2 ~ "TUTR",
        
        # ───────── otherwise keep existing label ───────────────────────────
        TRUE ~ UTR_type
      )
    ) %>% 
    ungroup() %>% 
    
    ## ───────── 2.  MIXED labelling per ensembl_ID ──────────
    group_by(ensembl_ID) %>% 
    mutate(
      UTR_type = case_when(
        any(UTR_type == "TUTR") ~ "MIXED",   # TRUE branch recycled automatically
        TRUE                     ~ UTR_type  # fallback
      )
    ) %>% 
    ungroup()
})

#split into iAPA and MIXED
MIXED_sites_list = lapply(iAPA_sites_list, function(df){
  df = df %>% filter(UTR_type == 'MIXED')
})

iAPA_sites_list = lapply(iAPA_sites_list, function(df){
  df = df %>% filter(UTR_type == 'iAPA')
})

# Process remaining 'UNKNOWN' sites
remaining_unknown_sites_list <- lapply(remaining_unknown_sites_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'UNKNOWN')
})

# Classify as 'ALE' or 'MIXED' based on criteria
sites_with_ALE_UTR_type_list <- lapply(remaining_unknown_sites_list, function(df) {
  df %>%
    group_by(ensembl_ID, splice_start) %>%
    mutate(num_unique_starts = n_distinct(start)) %>%
    ungroup() %>%
    group_by(ensembl_ID) %>%
    mutate(
      has_multiple_starts = any(num_unique_starts > 1),
      UTR_type = ifelse(has_multiple_starts, "MIXED", "ALE")
    ) %>%
    dplyr::select(-has_multiple_starts, -num_unique_starts) %>%
    ungroup()
})

# Separate 'ALE' sites
ALE_sites_list <- lapply(sites_with_ALE_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'ALE')
})

# Separate 'MIXED' sites
more_MIXED_sites_list <- lapply(sites_with_ALE_UTR_type_list, function(df) {
  df %>% dplyr::filter(UTR_type == 'MIXED')
})

# Combine 'MIXED' sites
MIXED_sites_list <- Map(rbind, MIXED_sites_list, more_MIXED_sites_list)

# Combine all sites into one list
all_sites_with_UTR_type <- Map(rbind, TUTR_sites_list, iAPA_sites_list)
all_sites_with_UTR_type <- Map(rbind, all_sites_with_UTR_type, ALE_sites_list)
all_sites_with_UTR_type <- Map(rbind, all_sites_with_UTR_type, MIXED_sites_list)

#replace the gene_name column with values from my annotation gene_name column (extracted from further annotation)
all_sites_with_UTR_type = lapply(all_sites_with_UTR_type,function(df){         
  df = df %>% mutate(
    gene_name = str_extract(
      further_annotation,
      "(?<=gene_name\\s)[^;]+"  # look-behind, then “not-;” characters
    )
  )
})

# Keep relevant columns and remove duplicate rows
all_sites_with_UTR_type <- lapply(all_sites_with_UTR_type, function(df) {
  df <- df[, c(2,1,3:17, 35, 43), drop = FALSE] %>% distinct()
  return(df)
})

#add new UTR_to_internal_shift value for the iAPA genes
all_sites_with_UTR_type_with_internal_to_UTR_shift = lapply(all_sites_with_UTR_type, function(df) {
  df <- df %>% 
    group_by(gene_name) %>% 
    mutate(
      internal_to_UTR_shift = {
        if (any(UTR_type == "iAPA")) {
          val_NO  <- change_in_usage[near_UTR_boundary == "NO"][1]
          val_YES <- change_in_usage[near_UTR_boundary == "YES"][1]
          if (length(val_NO) == 0 || length(val_YES) == 0) NA_real_
          else val_YES - val_NO
        } else {
          NA_real_          # numeric NA for non‑iAPA genes
        }
      }
    ) %>% 
    ungroup() %>% select(-near_UTR_boundary)
})

# ----------------------------
# Filter Significant and Control Genes for Barcharts
# ----------------------------

# Filter to isolate significant genes (genes with at least two significant PAS)
all_sites_with_UTR_type_with_internal_to_UTR_shift <- lapply(all_sites_with_UTR_type_with_internal_to_UTR_shift, function(df) {
  df <- df %>%
  mutate(is_significant = ifelse(!is.na(twostep_feature_padj) & 
                                   twostep_feature_padj < 0.05 & 
                                   abs(change_in_usage) > 0.1, TRUE, FALSE))

  # Count the number of significant PolyA sites per gene
  df <- df %>%
  group_by(gene_name) %>%
  mutate(sig_PAS_count = sum(is_significant, na.rm = TRUE)) %>%
  ungroup()
return(df)
})

# For genes with at least two significant PAS
sig_sites_with_UTR_type = lapply(all_sites_with_UTR_type_with_internal_to_UTR_shift,function(df){
  df_sig <- df %>%
  dplyr::filter(sig_PAS_count >= 2 & is_significant) %>%  # Keep only significant PAS
  group_by(gene_name) %>%
  arrange(desc(abs(change_in_usage))) %>%  # Order by absolute change in usage
  slice_head(n = 2) %>%  # Select top 2 PAS per gene
  ungroup() %>% dplyr::select(-is_significant)

  return(df_sig)
})


# Filter to isolate control genes (genes with stable PAS usage)
ctrl_sites_with_UTR_type <- lapply(all_sites_with_UTR_type_with_internal_to_UTR_shift, function(df) {
  df %>%
    dplyr::filter(mean_usage > 0.075 & abs(change_in_usage) < 0.075) %>%
    group_by(gene_name) %>%
    dplyr::filter(n() >= 2) %>%
    ungroup() %>% dplyr::select(-is_significant)
})

# ----------------------------
# Counting Genes per UTR Type and Plotting
# ----------------------------

# Function to count genes for each 'UTR_type'
count_genes_by_UTR_type <- function(df_list, target_UTR_type) {
  lapply(df_list, function(df) {
    df %>%
      group_by(gene_name) %>%
      dplyr::filter(any(UTR_type == target_UTR_type)) %>%
      ungroup() %>%
      summarise(num_genes = n_distinct(gene_name))
  })
}

# Function to safely extract 'num_genes' from each dataframe
extract_num_genes <- function(df) {
  if (nrow(df) == 0) {
    return(0)
  } else {
    return(df$num_genes)
  }
}

# Count number of genes belonging to each 'UTR_type' for significant genes
sig_TUTR_genes_count <- count_genes_by_UTR_type(sig_sites_with_UTR_type, "TUTR")
sig_ALE_genes_count <- count_genes_by_UTR_type(sig_sites_with_UTR_type, "ALE")
sig_MIXED_genes_count <- count_genes_by_UTR_type(sig_sites_with_UTR_type, "MIXED")
sig_iAPA_genes_count <- count_genes_by_UTR_type(sig_sites_with_UTR_type, "iAPA")

# Create a data frame with counts for each 'UTR_type' and each condition
UTR_gene_tally_list <- lapply(1:length(sig_TUTR_genes_count), function(i) {
  data.frame(
    UTR_type = c('TUTR', 'ALE', 'MIXED', 'iAPA'),
    total = c(
      extract_num_genes(sig_TUTR_genes_count[[i]]),
      extract_num_genes(sig_ALE_genes_count[[i]]),
      extract_num_genes(sig_MIXED_genes_count[[i]]),
      extract_num_genes(sig_iAPA_genes_count[[i]])
    ),
    condition = conditions[i]
  )
})

# Combine all tallies into one data frame for plotting
combined_tallies <- do.call(rbind, UTR_gene_tally_list)

# Convert 'UTR_type' to factor with specified levels for consistent plotting
combined_tallies$UTR_type <- factor(combined_tallies$UTR_type, levels = c('TUTR', 'ALE', 'MIXED', 'iAPA'))

# Define custom labels for the strip titles
strip_labels <- c(
  "ADMAi" = "ADMAi",
  "SDMAi" = "SDMAi",
  "SDMAi_ADMAi" = "SDMAi+ADMAi"
)

# Plot the bar chart for significant genes
sig_UTR_class_barchart <- ggplot(combined_tallies, aes(x = UTR_type, y = total, fill = UTR_type)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.7) +
  facet_wrap(~condition, scales = "free_x", labeller = labeller(condition = strip_labels)) +
  labs(
    x = 'APA Type',
    y = 'Number of Genes'
  ) +
  theme_classic() +
  scale_fill_brewer(palette = 'Set2') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgrey", color = "black", size = 1)
  )

# Save the bar chart plot to a file
ggsave(filename = file.path(output_dir, "sig_UTR_class_barchart_new.jpg"), 
       plot = sig_UTR_class_barchart, width = 10, height = 8, dpi = 300)

# ----------------------------
# Repeat Counting and Plotting for Control Genes
# ----------------------------

# Count number of genes belonging to each 'UTR_type' for control genes
ctrl_TUTR_genes_count <- count_genes_by_UTR_type(ctrl_sites_with_UTR_type, "TUTR")
ctrl_ALE_genes_count <- count_genes_by_UTR_type(ctrl_sites_with_UTR_type, "ALE")
ctrl_MIXED_genes_count <- count_genes_by_UTR_type(ctrl_sites_with_UTR_type, "MIXED")
ctrl_iAPA_genes_count <- count_genes_by_UTR_type(ctrl_sites_with_UTR_type, "iAPA")

# Create a data frame with counts for each 'UTR_type' and each condition
UTR_gene_tally_ctrl_list <- lapply(1:length(ctrl_TUTR_genes_count), function(i) {
  data.frame(
    UTR_type = c('TUTR', 'ALE', 'MIXED', 'iAPA'),
    total = c(
      extract_num_genes(ctrl_TUTR_genes_count[[i]]),
      extract_num_genes(ctrl_ALE_genes_count[[i]]),
      extract_num_genes(ctrl_MIXED_genes_count[[i]]),
      extract_num_genes(ctrl_iAPA_genes_count[[i]])
    ),
    condition = conditions[i]
  )
})

# Combine all tallies into one data frame for plotting
combined_tallies_ctrl <- do.call(rbind, UTR_gene_tally_ctrl_list)

# Convert 'UTR_type' to factor with specified levels for consistent plotting
combined_tallies_ctrl$UTR_type <- factor(combined_tallies_ctrl$UTR_type, levels = c('TUTR', 'ALE', 'MIXED', 'iAPA'))

# Plot the bar chart for control genes
ctrl_UTR_class_barchart <- ggplot(combined_tallies_ctrl, aes(x = UTR_type, y = total, fill = UTR_type)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.7) +
  facet_wrap(~condition, scales = "free_x", labeller = labeller(condition = strip_labels)) +
  labs(
    x = 'UTR Type',
    y = 'Number of Genes',
  ) +
  theme_classic() +
  scale_fill_brewer(palette = 'Set2') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgrey", color = "black", size = 1)
  )

# Save the bar chart plot to a file
ggsave(filename = file.path(output_dir, "ctrl_UTR_class_barchart.jpg"), 
       plot = ctrl_UTR_class_barchart, width = 10, height = 8, dpi = 300)

# ----------------------------
# Jitter Plot for Proximal to Distal Usage Shift and write CSVs to file
# ----------------------------

# Define functions to identify positive and negative strands
identify_positive_strands <- function(df) {
  df %>% dplyr::filter(strand == "+")
}

identify_negative_strands <- function(df) {
  df %>% dplyr::filter(strand == "-")
}

# Define functions to identify distal and proximal PolyA sites
identify_distal_pAs <- function(pos_df, neg_df) {
  pos_distal <- pos_df %>% group_by(gene_name) %>% dplyr::slice_max(order_by = start, n = 1, with_ties = FALSE)
  neg_distal <- neg_df %>% group_by(gene_name) %>% dplyr::slice_min(order_by = start, n = 1, with_ties = FALSE)
  bind_rows(pos_distal, neg_distal)
}

identify_proximal_pAs <- function(pos_df, neg_df) {
  pos_proximal <- pos_df %>% group_by(gene_name) %>% dplyr::slice_min(order_by = start, n = 1, with_ties = FALSE)
  neg_proximal <- neg_df %>% group_by(gene_name) %>% dplyr::slice_max(order_by = start, n = 1, with_ties = FALSE)
  bind_rows(pos_proximal, neg_proximal)
}

# **Step 1:** Identify positive and negative strand transcripts for all genes
pos_sites_with_UTR_type_anno_list <- lapply(all_sites_with_UTR_type_with_internal_to_UTR_shift, identify_positive_strands)
neg_sites_with_UTR_type_anno_list <- lapply(all_sites_with_UTR_type_with_internal_to_UTR_shift, identify_negative_strands)

# **Step 2:** Identify distal and proximal PolyA sites
anno_distal_pAs_list <- Map(identify_distal_pAs, pos_sites_with_UTR_type_anno_list, neg_sites_with_UTR_type_anno_list)
anno_proximal_pAs_list <- Map(identify_proximal_pAs, pos_sites_with_UTR_type_anno_list, neg_sites_with_UTR_type_anno_list)

# **Step 3:** Combine distal and proximal PolyA sites with 'pA_type' column
all_sites_with_UTR_type_anno_with_position_list <- Map(function(dist_df, prox_df) {
  dist_df$pA_type <- "distal pA"
  prox_df$pA_type <- "proximal pA"
  bind_rows(dist_df, prox_df)
}, anno_distal_pAs_list, anno_proximal_pAs_list)

# **Step 4:** Add 'condition' and information to each dataframe and write CSV files 
for (i in seq_along(all_sites_with_UTR_type_anno_with_position_list)) {
  all_sites_with_UTR_type_anno_with_position_list[[i]]$condition <- conditions[i]
}

setwd('/Users/llywelyngriffith/Documents/AZ_postdoc/quantseq_96hrs/Nobby_APA_analysis/dedup/github/with_intronic/TUTR_vs_ALE/CSV')

for (df in all_sites_with_UTR_type_anno_with_position_list) {
  # Extract the condition value from the 'condition' column
  condition_value <- unique(df$condition)
  
  # Construct the filename
  filename <- paste0("all_", condition_value, "_sites_with_UTR_type_info.csv")
  
  # Full file path
  filepath <- file.path(anno_bed_dir,'CSV', filename)
  
  # Write the data frame to CSV
  write_csv(df, filepath)
}

# ----------------------------
# Filter Significant and Control Genes for CSV files and Jitter plot
# ----------------------------

# Filter to isolate significant genes (genes with at least two significant PAS)
sig_sites_with_UTR_type_anno_with_position_list <- lapply(all_sites_with_UTR_type_anno_with_position_list, function(df) {
  df %>% dplyr::filter(sig_PAS_count >= 2)
})

# Filter to isolate control genes (genes with stable PAS usage)
ctrl_sites_with_UTR_type_anno_with_position_list <- lapply(all_sites_with_UTR_type_anno_with_position_list, function(df) {
  df %>%
    dplyr::filter(mean_usage > 0.075 & abs(change_in_usage) < 0.075) %>%
    group_by(gene_name) %>%
    dplyr::filter(n() >= 2) %>%
    ungroup()
})

for (df in sig_sites_with_UTR_type_anno_with_position_list) {
  # Extract the condition value from the 'condition' column
  condition_value <- unique(df$condition)
  
  # Construct the filename
  filename <- paste0("sig_", condition_value, "_sites_with_positional_info.csv")
  
  # Full file path
  filepath <- file.path(polyA_csv_dir,'positional_info', filename)
  
  # Write the data frame to CSV
  write_csv(df, filepath)
}

# **Step 5:** Extract specific UTR types for significant genes 
sig_TUTR_genes_list <- lapply(sig_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "TUTR")
  df_filtered$UTR_type <- "TUTR"  # Ensure consistency
  return(df_filtered)
})

sig_ALE_genes_list <- lapply(sig_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "ALE")
  df_filtered$UTR_type <- "ALE"
  return(df_filtered)
})

sig_MIXED_genes_list <- lapply(sig_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "MIXED")
  df_filtered$UTR_type <- "MIXED"
  return(df_filtered)
})

sig_iAPA_genes_list <- lapply(sig_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "iAPA")
  df_filtered$UTR_type <- "iAPA"
  return(df_filtered)
})

# **Step 6:** Combine all UTR types into one data frame per condition
sig_prox_dist_shift_list <- list()

for (i in seq_along(sig_sites_with_UTR_type_anno_with_position_list)) {
  # Extract data frames
  all_df <- sig_sites_with_UTR_type_anno_with_position_list[[i]]  # all UTR types
  TUTR_df <- sig_TUTR_genes_list[[i]]                             # TUTR
  ALE_df <- sig_ALE_genes_list[[i]]                               # ALE
  MIXED_df <- sig_MIXED_genes_list[[i]]                           # MIXED
  iAPA_df <- sig_iAPA_genes_list[[i]]                             # iAPA
  
  # Combine all UTR types into one data frame
  combined_df <- bind_rows(all_df, TUTR_df, ALE_df, MIXED_df, iAPA_df)
  
  # Identify distal and proximal PolyA sites
  pos_df <- combined_df %>% dplyr::filter(strand == "+")
  neg_df <- combined_df %>% dplyr::filter(strand == "-")
  
  distal_pAs <- identify_distal_pAs(pos_df, neg_df)
  proximal_pAs <- identify_proximal_pAs(pos_df, neg_df)
  
  # Add 'pA_type' and combine
  distal_pAs$pA_type <- "distal pA"
  proximal_pAs$pA_type <- "proximal pA"
  prox_dist_combined <- bind_rows(distal_pAs, proximal_pAs)
  
  # Calculate 'prox_to_dist_shift' per gene, UTR_type, and condition
  prox_dist_shift <- prox_dist_combined %>%
    group_by(gene_name, UTR_type, condition) %>%
    mutate(prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] - change_in_usage[pA_type == "proximal pA"]) %>%
    ungroup()
  
  # Remove outlier genes based on 'prox_to_dist_shift'
  prox_dist_shift_filtered <- prox_dist_shift %>% 
    dplyr::filter(abs(prox_to_dist_shift) > 0.2)
  
  # Append to the list
  sig_prox_dist_shift_list[[i]] <- prox_dist_shift_filtered
}

# **Step 7:** Combine all processed data into one data frame
combined_prox_dist_shift <- bind_rows(sig_prox_dist_shift_list)

#remove iAPA for prox to shift boxplot
combined_prox_dist_shift_no_iAPA = combined_prox_dist_shift %>% filter(UTR_type != 'iAPA')

# Duplicate the data frame and change UTR_type to 'all' in the duplicate
all_utr_df <- combined_prox_dist_shift_no_iAPA %>% 
  mutate(UTR_type = 'all')

# Combine the original data frame with the modified duplicate
combined_prox_dist_shift_no_iAPA <- bind_rows(combined_prox_dist_shift_no_iAPA, all_utr_df)

#filter to only keep All, TUTR, ALE and MIXED
combined_prox_dist_shift_no_iAPA = combined_prox_dist_shift_no_iAPA %>% filter(UTR_type == 'all' | UTR_type == 'TUTR' | UTR_type == 'ALE' | UTR_type == 'MIXED')

# **Step 8:** Ensure 'UTR_type' and 'pA_type' are factors with specified levels
combined_prox_dist_shift_no_iAPA$UTR_type <- factor(combined_prox_dist_shift_no_iAPA$UTR_type, levels = c('all','TUTR', 'ALE', 'MIXED'))

# Get the first four colors from the Set2 palette
set2_colors <- brewer.pal(8, "Set2")[1:3]

# Construct the color vector
color_vector <- c(
  '#ED6262',          # First color
  set2_colors,        # First four elements from Set2
  '#F3ADAD',          # Middle color
  set2_colors,        # Repeat first four elements from Set2
  '#9C2A2A',          # Last color
  set2_colors         # Repeat first four elements from Set2 again
)

# Create a unique identifier for each combination of UTR_type and condition
combined_prox_dist_shift_no_iAPA$group_id <- with(combined_prox_dist_shift_no_iAPA, interaction(UTR_type, condition))
# Check how many unique combinations there are
unique_combinations <- unique(combined_prox_dist_shift_no_iAPA$group_id)
# Create a mapping from combinations to colors
color_mapping <- setNames(color_vector, unique_combinations)
# Map the colors to your dataframe
combined_prox_dist_shift_no_iAPA$colors <- color_mapping[combined_prox_dist_shift_no_iAPA$group_id]

#add condition labels
no_iAPA_condition_labels <- c(
  "ADMAi"      = "ADMAi",
  "SDMAi"      = "SDMAi",
  "SDMAi_ADMAi" = "SDMAi+ADMAi"
)

jitter_plot = ggplot(combined_prox_dist_shift_no_iAPA, aes(x = UTR_type, y = prox_to_dist_shift, fill = colors)) +
  geom_violin(trim = FALSE, alpha = 0.4) + # Add violin plot to show density
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.3), 
              size = 1.5, aes(color = colors), alpha = 0.6) + # Add jitter points for visibility
  stat_summary(fun = median, geom = "point", 
               position = position_dodge(width = 0.3), 
               size = 6, shape = 22, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.8) +
  facet_wrap(~condition, scales = "free_x", labeller = labeller(condition = strip_labels)) +
  scale_fill_identity() +
  scale_color_identity() +
  labs(
    x = 'UTR Type',
    y = "Proximal to Distal Usage Shift",
    color = 'UTR Type',
    fill = 'UTR Type'
  ) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgrey", color = "black", size = 1)
  ) +
  ylim(-1.9, 1.9)


#separate into iAPA and non-iAPA
combined_prox_dist_shift_just_iAPA = combined_prox_dist_shift %>% filter(UTR_type == 'iAPA')

#create the named vector you pass to labeller()
condition_labels <- c(
  ADMAi = "ADMAi",
  SDMAi = "SDMAi",      
  SDMAi_ADMAi   = "SDMAi +\nADMAi"
)

# the 4-th colour of the 8-colour Set2 palette
set2_4 <- brewer.pal(8, "Set2")[4]      # "#E78AC3"
ggplot(combined_prox_dist_shift_just_iAPA,
       aes(x = UTR_type, y = internal_to_UTR_shift)) +
  
  geom_violin(trim  = FALSE,
              fill  = set2_4,   # pink fill
              colour = "black", # ← black border
              alpha = 0.4) +
  
  geom_jitter(width  = 0.20,
              size   = 1.5,
              colour = set2_4,
              alpha  = 0.6) +
  
  stat_summary(fun  = median, geom = "point",
               size = 6, shape = 22, fill = "white") +
  
  geom_hline(yintercept = 0, linetype = "dotted",
             colour = "black", size = 0.8) +
  
  facet_wrap(~condition,
             scales   = "free_x",
             labeller = labeller(condition = condition_labels)) +
  
  labs(x = "APA Type",
       y = "Internal to 3' UTR Shift") +
  
  theme_classic() +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    strip.text = element_text(size = 16),
    panel.border     = element_rect(colour = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "lightgrey",
                                    colour = "black", size = 1)
  ) +
  ylim(-1.75, 1.75)



# Save the jitter plot to a file
ggsave(filename = file.path(output_dir, "proximal_to_distal_usage_shift_TUTR_vs_ALE_jitter_plot.jpg"), 
       plot = jitter_plot, width = 10, height = 8, dpi = 300)

# ----------------------------
# Write Significant and Control Data Frames to CSV Files
# ----------------------------

# Function to write data frames to CSV files
write_df_list_to_files <- function(df_list, prefix) {
  # Categories corresponding to each dataframe in the list
  categories <- c("all", "TUTR", "ALE", "MIXED", "iAPA")
  
  # Check if the list has exactly five dataframes
  if(length(df_list) != length(categories)) {
    stop(paste("Expected", length(categories), "dataframes in the list, but got", length(df_list)))
  }
  
  # Iterate through the list and write each dataframe to a CSV file
  for(i in seq_along(df_list)) {
    # Construct the filename based on the category and condition
    filename <- paste0(prefix, "_", categories[i], ".csv")
    
    # Write the dataframe to a CSV file without row names
    write.csv(df_list[[i]], file = file.path(output_dir, filename), row.names = FALSE)
    
    # Optional: Print a message indicating successful write
    message(paste("Written:", filename))
  }
}

# Prepare significant data frames for writing
for (i in seq_along(sig_sites_with_UTR_type)) {
  condition <- conditions[i]
  
  # Create lists of data frames for each condition
  sig_list <- list(
    sig_sites_with_UTR_type[[i]],  # all
    sig_TUTR_genes_list[[i]],      # TUTR
    sig_ALE_genes_list[[i]],       # ALE
    sig_MIXED_genes_list[[i]],     # MIXED
    sig_iAPA_genes_list[[i]]       # iAPA
  )
  
  # Write to CSV files
  write_df_list_to_files(sig_list, paste0("sig_", condition))
}

#same for controls
ctrl_TUTR_genes_list <- lapply(ctrl_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "TUTR")
  df_filtered$UTR_type <- "TUTR"  # Ensure consistency
  return(df_filtered)
})

ctrl_ALE_genes_list <- lapply(ctrl_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "ALE")
  df_filtered$UTR_type <- "ALE"
  return(df_filtered)
})

ctrl_MIXED_genes_list <- lapply(ctrl_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "MIXED")
  df_filtered$UTR_type <- "MIXED"
  return(df_filtered)
})

ctrl_iAPA_genes_list <- lapply(ctrl_sites_with_UTR_type_anno_with_position_list, function(df) {
  df_filtered <- df %>% dplyr::filter(UTR_type == "iAPA")
  df_filtered$UTR_type <- "iAPA"
  return(df_filtered)
})

# Prepare control data frames for writing
for (i in seq_along(ctrl_sites_with_UTR_type)) {
  condition <- conditions[i]
  
  # Create lists of data frames for each condition
  ctrl_list <- list(
    ctrl_sites_with_UTR_type[[i]],  # all
    ctrl_TUTR_genes_list[[i]],      # TUTR
    ctrl_ALE_genes_list[[i]],       # ALE
    ctrl_MIXED_genes_list[[i]],     # MIXED
    ctrl_iAPA_genes_list[[i]]       # iAPA
  )
  
  # Write to CSV files
  write_df_list_to_files(ctrl_list, paste0("ctrl_", condition))
}

##top2 positional boxplot without iAPA
#first add prox to dist shift
all_sites_with_UTR_type_anno_with_prox_to_dist_shift = lapply(all_sites_with_UTR_type_anno_with_position_list,function(df) {
  # Identify distal and proximal PolyA sites
  pos_df <- df %>% dplyr::filter(strand == "+")
  neg_df <- df %>% dplyr::filter(strand == "-")
  
  distal_pAs <- identify_distal_pAs(pos_df, neg_df)
  proximal_pAs <- identify_proximal_pAs(pos_df, neg_df)
  
  # Add 'pA_type' and combine
  distal_pAs$pA_type <- "distal pA"
  proximal_pAs$pA_type <- "proximal pA"
  prox_dist_combined <- bind_rows(distal_pAs, proximal_pAs)
  
  # Calculate 'prox_to_dist_shift' per gene, UTR_type, and condition
  prox_dist_shift <- prox_dist_combined %>%
    group_by(gene_name, UTR_type, condition) %>%
    mutate(prox_to_dist_shift = change_in_usage[pA_type == "distal pA"] - change_in_usage[pA_type == "proximal pA"]) %>%
    ungroup()
  
  return(prox_dist_shift)
})

#remove iAPA from tables
all_sites_with_UTR_type_anno_with_prox_to_dist_shift = lapply(all_sites_with_UTR_type_anno_with_prox_to_dist_shift,function(df){
  df = df %>% filter(UTR_type != 'iAPA')
  return(df)
})

# Combine the shifted PAS data across conditions for plotting
combined_shifted_pas <- bind_rows(all_sites_with_UTR_type_anno_with_prox_to_dist_shift) 

# Set factor levels for plotting
combined_shifted_pas$condition <- factor(combined_shifted_pas$condition, levels = c('ADMAi', 'SDMAi' ,'SDMAi_ADMAi'))
combined_shifted_pas$pA_type <- factor(combined_shifted_pas$pA_type, levels = c('proximal pA', 'distal pA'))

# Step 24: Plot the proximal to distal usage shift
# Remove top 5% highest absolute shifts in each group
filtered_data <- combined_shifted_pas %>%
  group_by(condition) %>%
  mutate(abs_shift = abs(prox_to_dist_shift)) %>%
  arrange(desc(abs_shift)) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(rank > 0.1 * max(rank)) %>%
  ungroup() %>%
  dplyr::select(-abs_shift, -rank) %>% # Remove temporary columns
  group_by(gene_name,condition) %>% dplyr::slice(1) %>% ungroup() # only need one row per gene

#define colour scheme
top2_colours = c('#ED6262','#F3ADAD','#9C2A2A')

# Create the boxplot
shift_boxplot <- ggplot(filtered_data, aes(x = condition, y = prox_to_dist_shift, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = top2_colours,labels = c("ADMAi", "SDMAi", "SDMAi+ADMAi")) +
  ylab('Proximal to Distal Usage Shift') +
  scale_x_discrete(labels=c("ADMAi" = "ADMAi", "SDMAi" = "SDMAi","SDMAi_ADMAi" = "SDMAi+ADMAi")) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.8) +
  theme_classic() +
  ylim(-0.48, 0.48) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 20)
  ) +
  labs(fill = "Condition")
