#!/bin/sh

#### this script generates fasta files from input bed files

input_dir="ulelab/PRMT-APA/Data/Figure_4_S7/Figure_4_S7/bed/APA_3seq/1000nt_window"  # change to boxplot_regions/ or /200nt_window to generate those respective fasta files
outdir_dir="ulelab/PRMT-APA/Data/Figure_4_S7/Figure_4_S7/fasta/APA_3seq/1000nt_window"

# Iterate over "*.bed" files within the input directory
for bedfile in "$input_dir"/*.bed; do
    # Extract the base filename without extension
    filename=$(basename "$bedfile" .bed)
    
    # Define the output file path
    outputfile="$output_dir/${filename}.fasta"
    
   # Run bedtools to convert bed to fasta
    bedtools getfasta -fi /path/to/GRCh38.primary_assembly.genome.fa -bed "$bedfile" -fo "$outputfile" -s -name
done


