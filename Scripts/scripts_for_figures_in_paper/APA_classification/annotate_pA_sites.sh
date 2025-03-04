#!/bin/sh

# Description:
# This script processes CSV files in a specified input folder to create BED files,
# then uses bedtools to annotate the BED files against a GTF annotation file. poly(A) sites that do not overlap with any gene within the gtf file are removed and stored in a 'no_annotation' bed file.
# Ensure that bedtools is installed and available in your system's PATH.

# ----------------- User Configuration -----------------

# Set the input folder containing the CSV files
# Change this to the path where your CSV files are located
input_folder="/path/to/your/csv/files/"
input_folder="/ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/CSV_files/Fig1
/DMAi_timecourse/"

# Set the path to the GTF annotation file
# Change this to the path where your GTF annotation file is located
gtf_annotation="/path/to/your/gencode.v45.annotation.gtf"

# Set the output folder where the results will be stored
# Change this to your desired output directory e.g:
output_folder="/ulelab/repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/bed_files/"

# ------------------------------------------------------

# Iterate through each CSV file in the input folder matching the pattern
for input_file in "$input_folder"*"_with"*"_change_in_usage.csv"; do
    # Extract the base name of the file (excluding path and specific suffix)
    base_name=$(basename "$input_file" "_with_change_in_usage.csv")
    
    echo "Processing file: $input_file"
    echo "Base name: $base_name"
    
    # Define the intermediate BED file path
    bed_file="${output_folder}${base_name}_coords.bed"
    
    # Convert the CSV file to BED format
    # Adjust the column numbers in awk if necessary to match your CSV file structure
    # The example assumes columns: 13 (chr), 14 (start), 15 (end), 1 (gene_name), 16 (strand)
    awk -F ',' 'NR>1 {print $13 "\t" $14 "\t" $15 "\t" "." "\t" $1 "\t" $16}' "$input_file" | sort -k1,1 -k2,2n > "$bed_file"
    
    echo "Created BED file: $bed_file"
    
    # Define the output files for the bedtools commands
    annotated_bed="${output_folder}${base_name}_coords_annotated.bed"
    no_annotation_bed="${output_folder}${base_name}_coords_no_annotation.bed"
    
    # Run bedtools window to annotate the BED file
    # Adjust the bedtools parameters as needed
    bedtools window -l 100 -r 1 -sw -a "$bed_file" -b "$gtf_annotation" -sm > "$annotated_bed"
    
    echo "Created annotated BED file: $annotated_bed"
    
    # Run bedtools window with -v option to get entries with no annotation
    bedtools window -l 100 -r 1 -sw -a "$bed_file" -b "$gtf_annotation" -sm -v > "$no_annotation_bed"
    
    echo "Created BED file with no annotation: $no_annotation_bed"
    
done

echo "Processing complete."
