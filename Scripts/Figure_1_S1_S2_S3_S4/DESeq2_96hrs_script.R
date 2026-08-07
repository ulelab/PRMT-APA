##############################################################################
# 0. Load Required Packages
##############################################################################
library(tidyverse)       # Data manipulation
library(tximport)        # Importing transcript-level estimates
library(DESeq2)          # Differential expression analysis
library(apeglm)          # Shrinking log fold changes
library(org.Hs.eg.db)    # Human gene annotations
library(GenomicFeatures) # Creating TxDb from GTF/GFF

##############################################################################
# 1. Read and Prepare Sample Information
##############################################################################

# Read sample metadata
tso_seq_sample_info <- read_csv('Data/Figure_1_S1_S2_S3_S4/DESeq2_96hrs_sample_info.csv') %>%
  as.data.frame()  # Convert to data frame for DESeq2

# Make 'Treatment' a factor and set 'DMSO_96hrs' as the reference level
tso_seq_sample_info$Treatment <- factor(tso_seq_sample_info$Treatment)
tso_seq_sample_info$Treatment <- relevel(tso_seq_sample_info$Treatment, ref = "DMSO_96hrs")

# Base path relative to the *current working directory*
base_counts_dir <- "Data/Figure_1_S1_S2_S3_S4/salmon_files_for_DESeq2/96hrs"

# Create vector of Salmon quant.sf file paths and assign sample names
tso_seq_sample_files <- file.path(
  base_counts_dir,
  paste0(tso_seq_sample_info$Treatment, tso_seq_sample_info$Replicate),
  "quant.sf"
)
names(tso_seq_sample_files) <- paste0(tso_seq_sample_info$Treatment, '_', tso_seq_sample_info$Replicate)

# Create named vector for sample files (optional, but nice for clarity)
sample_list <- c(
  'DMSO_96hrs_1','DMSO_96hrs_2','DMSO_96hrs_3',
  'SDMAi_96hrs_1','SDMAi_96hrs_2','SDMAi_96hrs_3',
  'ADMAi_96hrs_1','ADMAi_96hrs_2','ADMAi_96hrs_3',
  'SDMAi_ADMAi_96hrs_1','SDMAi_ADMAi_96hrs_2','SDMAi_ADMAi_96hrs_3'
)
names(tso_seq_sample_files) <- sample_list

##############################################################################
# 2. Create tx2gene Mapping
##############################################################################
# Create a TxDb object from GTF or GFF
txdb <- makeTxDbFromGFF(
  file='/Users/k2362866/Documents/AZ_postdoc/CFIM25_SAM68_ELAVL1_siRNA_3_seq/Nobby_APA_analysis/own_atlas/TUTR_vs_ALE/gencode.v45.annotation.gtf'
)

# Build tx2gene map by linking GENEID to TXNAME
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene_map <- df[, 2:1]
# Remove the version suffix from TXNAME (e.g. ENSG00001.1 -> ENSG00001)
tx2gene_map$TXNAME <- sub("\\..*", "", tx2gene_map$TXNAME)

##############################################################################
# 3. Import Salmon Counts with tximport
##############################################################################
tximport_data <- tximport(
  files            = tso_seq_sample_files,
  type             = 'salmon',
  tx2gene          = tx2gene_map,
  countsFromAbundance = "no",
  ignoreTx         = TRUE
)

# Round counts to integers (required by DESeq2)
tximport_data$counts <- round(tximport_data$counts)

##############################################################################
# 4. Create DESeqDataSet and Filter
##############################################################################
tso_seq_dds <- DESeqDataSetFromMatrix(
  countData = tximport_data$counts,
  colData   = tso_seq_sample_info,
  design    = ~ Treatment
)

# Filter out low-count genes
keep_genes <- rowSums(counts(tso_seq_dds)[, 1:3] >= 5) >= 3 |
  rowSums(counts(tso_seq_dds)[, 4:6] >= 5) >= 3 |
  rowSums(counts(tso_seq_dds)[, 7:9] >= 5) >= 3 |
  rowSums(counts(tso_seq_dds)[, 10:12] >= 5) >= 3

tso_seq_dds <- tso_seq_dds[keep_genes, ]

##############################################################################
# 5. DESeq Analysis
##############################################################################
tso_seq_dds <- DESeq(tso_seq_dds)

# Perform variance stabilizing transformation for PCA
vsd <- vst(tso_seq_dds, blind = FALSE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = 'Treatment', returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Define custom colors and labels for treatments
my_colors <- c("DMSO_96hrs" = '#8A8A8A', "SDMAi_96hrs" = '#F3ADAD',"ADMAi_96hrs" = '#ED6262', "SDMAi_ADMAi_96hrs" = '#9C2A2A')

custom_labels <- c("DMSO_96hrs" = "DMSO 96hrs", 
                   "SDMAi_96hrs" = "SDMAi 96hrs", 
                   "ADMAi_96hrs" = "ADMAi 96hrs", 
                   "SDMAi_ADMAi_96hrs" = "SDMAi + ADMAi 96hrs")

# Plot PCA using ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = my_colors, labels = custom_labels) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(size = 16))
