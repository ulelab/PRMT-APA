#!/bin/sh

#### this script generates fasta files from input bed files

# Iterate over directories within the "bed" directory
for dir in /ulelab/PRMT-APA/Scripts/scripts_for_figures_in_paper/bed_files/200nt_window/DMAi_cancer_line_panel/*; do
    # Extract the directory name
    dir_name=$(basename "$dir")
    
    # Iterate over "*.bed" files within the directory
    for bedfile in "$dir"/*.bed; do
        filename=$(basename "$bedfile" .bed)
        outputfile="fasta/200nt_window/DMAi_cancer_line_panel/$dir_name/${filename}.fasta" #example directory, change to where you want to save fasta files 
        bedtools getfasta -fi /path/to/genome_fasta_file/GRCh38.primary_assembly.genome.fa -bed "$bedfile" -fo "$outputfile" -s -name
    done
done
