#!/bin/sh

#### this script generates fasta files from input bed files

# Iterate over directories within the "bed" directory
for dir in /ulelab/PRMT-APA/Data/Figure_3_S5_S6/bed/APA_3seq/1000nt_window/*; do
    # Extract the directory name
    dir_name=$(basename "$dir")
    
    # Iterate over "*.bed" files within the directory
    for bedfile in "$dir"/*.bed; do
        filename=$(basename "$bedfile" .bed)
        outputfile="/ulelab/PRMT-APA/Data/Figure_3_S5_S6/fasta/1000nt_window/$dir_name/${filename}.fasta" 
        bedtools getfasta -fi /path/to/genome_fasta_file/GRCh38.primary_assembly.genome.fa -bed "$bedfile" -fo "$outputfile" -s -name
    done
done
