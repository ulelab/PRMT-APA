# Description:
# This script uses the GTFtools package to generate bed files containing locations of splice sites and UTR regions

python GTFtools/gtftools.py -q ../../Data/Figure_1_S1_S2_S3/classification_of_APA_bed_files/splice_regions.bed /path/to/gencode.v45.annotation.gtf

python GTFtools/gtftools.py -u ../../Data/Figure_1_S1_S2_S3/classification_of_APA_bed_files/UTR_regions_gencode_v45.bed /path/to/gencode.v45.annotation.gtf
