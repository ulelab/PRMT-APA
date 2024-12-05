#!/bin/sh

#SBATCH -c 8

# Description:
# This script uses the GTFtools package to generate bed files containing locations of splice sites and 3' UTR regions

python GTFtools/gtftools.py -q splice_regions.bed /camp/lab/ulej/home/users/griffil2/RNA_seq/martina/TUTR_vs_ALE/gencode.v45.annotation.gtf

python GTFtools/gtftools.py -u UTR_regions_gencode_v45.bed /camp/lab/ulej/home/users/griffil2/RNA_seq/martina/TUTR_vs_ALE/gencode.v45.annotation.gtf
