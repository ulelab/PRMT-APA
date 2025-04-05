library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(data.table)

#### 1) make the 150nt bed files
# Specify the directory
bed_dir <- "/Users/llywelyngriffith/Documents/AZ_postdoc/files_for_charlotte/github/motif_and_clip_analysis/DMAi_vs_non_DMAi/bed/"

# Get full file paths of all .bed files in that directory
bed_file_paths <- list.files(
  path = bed_dir,
  pattern = "sites\\.bed$",  # Corrected pattern
  full.names = TRUE
)

# Specify desired column names
bed_file_colnames <- c("chr", "start", "end", "gene_name", "feature_id", "strand")

# Read each .bed file into R, storing them in a list
bed_files <- lapply(bed_file_paths, function(file) {
  read.delim(
    file,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = bed_file_colnames
  )
})

#read in whole pA atlas
atlas_colnames = c('chr','start', 'end','feature_id','score','strand','gene_symbol','gene_name','region')

pA_atlas = read.delim("/Users/llywelyngriffith/Documents/AZ_postdoc/Shaun_cell_lines/Nobby_APA_analysis/common_atlas/resequenced/dedup/merged_polya.filteredunique.annotated.bed",header = FALSE,
                      stringsAsFactors = FALSE,
                      col.names = atlas_colnames)

#create copy
bed_files_with_original_coordinates = bed_files

#replace start and end window values with values of pA coordinates, taken from pA atlas
for(i in seq_along(bed_files_with_original_coordinates)) {
  
  bed_files_with_original_coordinates[[i]] <- bed_files_with_original_coordinates[[i]] %>%
    # Left join pA_atlas by feature_id and gene_name
    left_join(pA_atlas %>% 
                dplyr::select(start, end, feature_id, gene_name),
              by = c("feature_id", "gene_name")) %>%
    
    # Overwrite the 'start' and 'end' columns with the ones from pA_atlas
    mutate(
      start = start.y,
      end   = end.y
    ) %>%
    
    # Remove the temporary columns from the join and put columns in correct order
    dplyr::select(-start.x, -end.x, -start.y, -end.y) %>% dplyr::select(chr, start, end, gene_name, feature_id, strand)
}

#write original coordinate tables to file
# Get the base names (no path) of the input files
bed_file_names <- basename(bed_file_paths)           
bed_file_stems <- sub("\\.bed$", "", bed_file_names) 
output_dir = "/Users/llywelyngriffith/Documents/AZ_postdoc/files_for_charlotte/github/motif_and_clip_analysis/DMAi_vs_non_DMAi/bed/relaxed_thresholds/SRA/sites"

#    write each data frame to a new .bed file
for (i in seq_along(bed_files_with_original_coordinates)) {
  
  # Build a new filename 
  new_bed_name <- paste0(bed_file_stems[i], ".bed")
  
  # Choose the output path; here we save alongside the original directory
  new_bed_path <- file.path(output_dir, new_bed_name)
  
  # Write the data frame to a .bed file
  write.table(
    bed_files_with_original_coordinates[[i]],
    file = new_bed_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
}

# generate window
APA_sites_with_window_list = lapply(bed_files_with_original_coordinates, function(df) {
  df <- df %>%
    mutate(
      start = if_else(strand == "+", end - 150, start),
      end   = if_else(strand == "+", end, start + 150)
    )
  return(df)
})

#output dir
output_dir = "/Users/llywelyngriffith/Documents/AZ_postdoc/files_for_charlotte/github/motif_and_clip_analysis/DMAi_vs_non_DMAi/bed/relaxed_thresholds/SRA/150_upstream"

# After generating 'APA_sites_with_window_list',
#    write each data frame to a new .bed file
for (i in seq_along(APA_sites_with_window_list)) {
  
  # Build a new filename by appending "_150upstream"
  new_bed_name <- paste0(bed_file_stems[i], "_150upstream.bed")
  
  # Choose the output path; here we save alongside the original directory
  new_bed_path <- file.path(output_dir, new_bed_name)
  
  # Write the data frame to a .bed file
  write.table(
    APA_sites_with_window_list[[i]],
    file = new_bed_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

}

###### 2) make proximal to distal region bed files and whole gene bed files
# create bed files of proximal to distal region
control_sites = rbind(bed_files_with_original_coordinates[[1]],bed_files_with_original_coordinates[[2]])
DMAi_sites = rbind(bed_files_with_original_coordinates[[3]],bed_files_with_original_coordinates[[4]])
JTE607_control_sites = rbind(bed_files_with_original_coordinates[[5]],bed_files_with_original_coordinates[[6]])
JTE607_sites = rbind(bed_files_with_original_coordinates[[7]],bed_files_with_original_coordinates[[8]])
non_DMAi_sites = rbind(bed_files_with_original_coordinates[[9]],bed_files_with_original_coordinates[[10]])
siCFIM25_sites = rbind(bed_files_with_original_coordinates[[11]],bed_files_with_original_coordinates[[12]])

#turn back into list
all_sites = list(control_sites,DMAi_sites,JTE607_control_sites,JTE607_sites,non_DMAi_sites,siCFIM25_sites) 
proximal_sites = list(bed_files_with_original_coordinates[[2]],bed_files_with_original_coordinates[[4]],bed_files_with_original_coordinates[[6]],bed_files_with_original_coordinates[[8]],bed_files_with_original_coordinates[[10]],bed_files_with_original_coordinates[[12]]) 

#read in gencode file
gencode_gtf_path='/Users/llywelyngriffith/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/own_atlas/TUTR_vs_ALE/gencode.v45.annotation.gtf'

# Import the GTF file
gtf <- import(gencode_gtf_path)

# Filter for exon features
exons <- gtf[gtf$type == "exon"]

# Convert to a data frame
exons_df <- as.data.frame(exons)

# Clean up columns: trim whitespace and ensure correct types
exons_df <- exons_df %>%
  mutate(
    gene_name = trimws(gene_name),
    strand = trimws(strand),
    seqnames = as.character(seqnames)
  )

# Extract the first exon per gene_name and strand
first_exons <- exons_df %>%
  group_by(gene_name, strand) %>%
  summarize(
    chr = dplyr::first(seqnames),  # Assuming all exons for a gene are on the same chromosome
    first_exon_coord = if_else(
      strand == "+",
      min(start),  # For '+' strand, use the smallest start
      max(end)     # For '-' strand, use the largest end
    ) ,
    .groups = "drop"
  )

#only keep one row per gene_name/strand combo
first_exons = first_exons %>% distinct(gene_name, strand, .keep_all = TRUE)

# Create a GRanges object for all exons
exons_gr <- GRanges(
  seqnames = exons_df$seqnames,
  ranges = IRanges(start = exons_df$start, end = exons_df$end),
  gene_name = exons_df$gene_name,
  strand = exons_df$strand
)

# Define a helper function to calculate exons_length
calculate_exons_length <- function(gene, strand, chr, region_start, region_end, exons_gr) {
  # 1) Subset exons_gr to just this gene (and correct strand)
  sub_exons_gr <- exons_gr[exons_gr$gene_name == gene & strand(exons_gr) == strand]
  if (length(sub_exons_gr) == 0) {
    return(0)
  }
  
  # 2) Now define the region
  region_gr <- GRanges(
    seqnames = chr,
    ranges   = IRanges(start = region_start, end = region_end),
    strand   = strand
  )
  
  # 3) Overlap *only* that gene’s exons
  overlapping_exons <- subsetByOverlaps(sub_exons_gr, region_gr, ignore.strand = FALSE)
  if (length(overlapping_exons) == 0) {
    return(0)
  }
  
  # 4) Clip and merge
  clipped_exons <- pintersect(overlapping_exons, region_gr, ignore.strand = FALSE)
  merged_exons  <- reduce(clipped_exons, ignore.strand = FALSE)
  
  sum(width(merged_exons))
}

## a) make the proximal to distal bed files
# Create a new list of dataframes
prox_to_dist_region_list <- lapply(all_sites, function(df) {
  
  # Summarize each gene_name (+ strand, + chr) to find the min/max coordinates
  # for the "distal region." The logic is the same as before except we also
  # keep track of chr and group by gene_name, strand, and chr so we can
  # pass them to calculate_exons_length().
  
  df2 <- df %>%
    group_by(gene_name, strand, chr) %>%
    summarize(
      UTR_region_start = if_else(
        strand == "+",
        min(end),    # for plus-strand, region starts at min(end)
        min(start)   # for minus-strand, region starts at min(start)
      ),
      UTR_region_end = if_else(
        strand == "+",
        max(end),    # for plus-strand, region ends at max(end)
        max(start)   # for minus-strand, region ends at max(start)
      ),
      .groups = "drop"
    ) %>% ungroup() %>% group_by(gene_name) %>% dplyr::slice(1) %>% ungroup()

  # Now compute the sum of merged exons within [UTR_region_start, UTR_region_end].
  # We'll use rowwise() so each row calls the 'calculate_exons_length' function 
  # with the appropriate region.
  
  df2 <- df2 %>%
    rowwise() %>%
    mutate(
      UTR_length = calculate_exons_length(
        gene   = gene_name,
        strand = strand,
        chr    = chr,
        region_start = UTR_region_start,
        region_end   = UTR_region_end,
        exons_gr     = exons_gr   # from your global environment
      )
    ) %>%
    ungroup()
  
  # Finally, select + reorder columns to match your desired .bed format
  df2 <- df2 %>%
    dplyr::select(chr, UTR_region_start, UTR_region_end, gene_name, UTR_length, strand) 
  
  return(df2)
})


output_dir = "/Users/llywelyngriffith/Documents/AZ_postdoc/files_for_charlotte/github/motif_and_clip_analysis/DMAi_vs_non_DMAi/bed/SRA/UTR_regions/test"

# After generating 'prox_to_dist_region_list',
#    write each data frame to a new .bed file
prefixes <- c(
  "control_",
  "DMAi_",
  "JTE607_control_",
  "JTE607_",
  "non_DMAi_",
  "siCFIM25_"
)

for (i in seq_along(prox_to_dist_region_list)) {
  
  # Build the new filename using the manually specified prefix
  new_bed_name <- paste0(prefixes[i], "sites_UTRregions.bed")
  
  # Full path to save the file
  new_bed_path <- file.path(output_dir, new_bed_name)
  
  # Write the data frame to a .bed file
  write.table(
    prox_to_dist_region_list[[i]],
    file = new_bed_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

## b) make the 1st exon to proximal site bed files
prox_to_exon_start_list <- lapply(proximal_sites, function(df) {

  # Join with first_exons to get first_exon_coord
  df_joined <- df %>%
    dplyr::left_join(first_exons %>% dplyr::select(gene_name, strand, first_exon_coord),
                     by = c("gene_name", "strand"))
  
  # Check for any NAs in first_exon_coord
  if(any(is.na(df_joined$first_exon_coord))) {
    warning("Some gene_name and strand combinations did not match in first_exons")
  }
  
  # Replace 'start' or 'end' based on strand
  df_joined <- df_joined %>%
    mutate(
      new_start = if_else(strand == "+", first_exon_coord, start),
      new_end = if_else(strand == "-", first_exon_coord, end)
    )
  
  # Calculate exons_length for each row
  df_joined <- df_joined %>%
    rowwise() %>%
    mutate(
      exons_length = if (strand == "+") {
        calculate_exons_length(
          gene = gene_name,
          strand = strand,
          chr = chr,
          region_start = first_exon_coord,
          region_end = end,
          exons_gr = exons_gr
        )
      } else if (strand == "-") {
        calculate_exons_length(
          gene = gene_name,
          strand = strand,
          chr = chr,
          region_start = start,
          region_end = first_exon_coord,
          exons_gr = exons_gr
        )
      } else {
        0  # Handle any unexpected strand values
      }
    ) %>%
    ungroup()
  
  # Replace 'start' or 'end' with the new values
  df_joined <- df_joined %>%
    mutate(
      start = if_else(strand == "+", new_start, start),
      end = if_else(strand == "-", new_end, end)
    )
  
  # Replace 'feature_id' with 'exons_length' 
  df_final <- df_joined %>%
    dplyr::select(chr, start, end, gene_name, exons_length, strand)  # Reorder columns
  
  return(df_final)
})

output_dir = "/Users/llywelyngriffith/Documents/AZ_postdoc/files_for_charlotte/github/motif_and_clip_analysis/DMAi_vs_non_DMAi/bed/SRA/first_exon_to_pPA/"

# 2) After generating 'prox_to_exon_start_list',
#    write each data frame to a new .bed file
prefixes <- c(
  "control_",
  "DMAi_",
  "JTE607_control_",
  "JTE607_",
  "non_DMAi_",
  "siCFIM25_"
)

for (i in seq_along(prox_to_exon_start_list)) {
  
  # Build the new filename using the manually specified prefix
  new_bed_name <- paste0(prefixes[i], "sites_first_exon_to_pPA.bed")
  
  # Full path to save the file
  new_bed_path <- file.path(output_dir, new_bed_name)
  
  # Write the data frame to a .bed file
  write.table(
    prox_to_exon_start_list[[i]],
    file = new_bed_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

