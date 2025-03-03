#!/bin/sh

# Iterate over directories within the "bed" directory
for dir in bed/new_dev_0_1/test/200/*; do
    # Extract the directory name
    dir_name=$(basename "$dir")
    
    # Iterate over "*.bed" files within the directory
    for bedfile in "$dir"/*.bed; do
        filename=$(basename "$bedfile" .bed)
        outputfile="fasta/new_dev_0_1/test/200/$dir_name/${filename}.fasta"
        bedtools getfasta -fi /camp/lab/ulej/home/users/griffil2/RNA_seq/quantseq_96hrs/Nobby_APA_analysis/genomes/GRCh38.primary_assembly.genome.fa -bed "$bedfile" -fo "$outputfile" -s -name
    done
done