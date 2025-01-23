# Description:
# This script uses the GTFtools package to generate bed files containing locations of splice sites and 3' UTR regions

python GTFtools/gtftools.py -q splice_regions.bed gencode.v45.annotation.gtf

python GTFtools/gtftools.py -u UTR_regions_gencode_v45.bed gencode.v45.annotation.gtf
