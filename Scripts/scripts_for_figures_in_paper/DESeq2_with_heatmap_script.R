# Load required packages
library(tidyverse)       # For data manipulation and visualization
library(tximport)        # For importing transcript-level estimates
library(DESeq2)          # For differential expression analysis
library(apeglm)          # For log fold change shrinkage
library(AnnotationDbi)   # For annotation functions
library(org.Hs.eg.db)    # For human gene annotations
library(GenomicFeatures) # For creating TxDb objects
library(ggplot2)         # For figures
library(RColorBrewer)    # For figures
library(pheatmap)        # For figures

#### DMAi timecourse DESeq2
# Set working directory - e.g DMAi timecourse
setwd('../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/')

# Read in sample metadata e.g
timecourse_sample_info <- read_csv('DESeq2_tables/deseq2_timecourse_sample_info.csv') %>%
  mutate(Treatment = factor(Treatment))

# Create vector of Salmon quant.sf file paths and assign sample names
timecourse_sample_files <- paste0(timecourse_sample_info$Treatment, timecourse_sample_info$Replicate, '/quant.sf')
names(timecourse_sample_files) <- paste0(timecourse_sample_info$Treatment, '_', timecourse_sample_info$Replicate)

# Create tx2gene mapping using GTF annotation
txdb = makeTxDbFromGFF(file='gencode.v45.annotation.gtf')
k = keys(txdb, keytype = "GENEID")
df = select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene_map = df[, 2:1]
tx2gene_map$TXNAME = sub("\\..*", "", tx2gene_map$TXNAME)

# Import transcript-level estimates and summarize to gene-level counts
txi <- tximport(files= timecourse_sample_files,
                type = 'salmon',
                tx2gene = tx2gene_map,
                countsFromAbundance="no",
                ignoreTx = TRUE)

# Round counts to integers (DESeq2 requires integer counts)
txi$counts <- round(txi$counts)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = txi$counts,
                              colData = timecourse_sample_info,
                              design = ~ Treatment)

# Filter out low-count genes (keep genes with counts >=5 in at least two samples)
keep <- rowSums(counts(dds) >= 5) >= 2
dds <- dds[keep, ]

# Run DESeq normalization and estimation
dds <- DESeq(dds)

# Perform variance stabilizing transformation for PCA
vsd <- vst(dds, blind = FALSE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = 'Treatment', returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Define custom colors and labels for treatments
my_colors <- c("DMSO_24hrs" = "#e4e4e4", "DMSO_48hrs" = "#c4c2c2", "DMSO_72hrs" = "#a7a6a6", "DMSO_96hrs" = "#606060",
               "SDMAi_ADMAi_24hrs" = "#f68f60", "SDMAi_ADMAi_48hrs" = "#ff2d2d", "SDMAi_ADMAi_72hrs" = "#da372a",
               "SDMAi_ADMAi_96hrs" = "#9C2A2A")
custom_labels <- c("DMSO_24hrs" = "DMSO 24hrs", "DMSO_48hrs" = "DMSO 48hrs", "DMSO_72hrs" = "DMSO 72hrs", "DMSO_96hrs" = "DMSO 96hrs",
                   "SDMAi_ADMAi_24hrs" = "DMAi 24hrs", "SDMAi_ADMAi_48hrs" = "DMAi 48hrs", "SDMAi_ADMAi_72hrs" = "DMAi 72hrs",
                   "SDMAi_ADMAi_96hrs" = "DMAi 96hrs")

# Plot PCA using ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colors, labels = custom_labels) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(size = 16))

# Differential expression analysis for each time point
time_points <- c("24hrs", "48hrs", "72hrs", "96hrs")

# Initialize list to store results
results_list <- list()

for (tp in time_points) {
  # Subset samples for the current time point
  treatments_tp <- c(paste0("DMSO_", tp), paste0("SDMAi_ADMAi_", tp))
  samples_tp <- timecourse_sample_info$Treatment %in% treatments_tp
  dds_tp <- dds[, samples_tp]
  
  # Relevel Treatment factor to set control as reference
  dds_tp$Treatment <- droplevels(dds_tp$Treatment)
  dds_tp$Treatment <- relevel(dds_tp$Treatment, ref = paste0("DMSO_", tp))
  
  # Identify samples for each condition within dds_tp
  condition_levels <- levels(dds_tp$Treatment)
  condition_A <- condition_levels[1]  # Reference level (e.g., "DMSO_24hrs")
  condition_B <- condition_levels[2]  # Treatment level (e.g., "SDMAi_ADMAi_24hrs")
  
  samples_condition_A <- colData(dds_tp)$Treatment == condition_A
  samples_condition_B <- colData(dds_tp)$Treatment == condition_B
  
  # Perform filtering: retain genes with counts ≥5 in at least two samples within the same condition
  keep_tp <- (rowSums(counts(dds_tp)[, samples_condition_A] >= 5) >= 2) |
    (rowSums(counts(dds_tp)[, samples_condition_B] >= 5) >= 2)
  dds_tp <- dds_tp[keep_tp, ]
  
  # Run DESeq
  dds_tp <- DESeq(dds_tp)
  
  # Get results with adjusted p-values
  res <- results(dds_tp, contrast = c("Treatment", condition_B, condition_A), alpha = 0.05)
  
  # Shrink LFC estimates using apeglm
  resLFC <- lfcShrink(dds_tp, coef = 2, type = "apeglm")
  
  # Convert results to data frame and add gene IDs
  resLFC_df <- as.data.frame(resLFC) %>%
    rownames_to_column(var = "ensgene")
  
  # Remove version numbers from Ensembl IDs
  resLFC_df$ensgene <- sub("\\..*", "", resLFC_df$ensgene)
  
  # Annotate genes with gene symbols
  gene_annotations <- select(org.Hs.eg.db, keys = resLFC_df$ensgene, columns = "SYMBOL", keytype = "ENSEMBL")
  colnames(gene_annotations) <- c("ensgene", "gene_symbol")
  resLFC_df_anno <- left_join(resLFC_df, gene_annotations, by = "ensgene")
  
  # Remove rows with NA gene_symbol and ensure uniqueness
  resLFC_df_anno <- resLFC_df_anno %>%
    filter(!is.na(gene_symbol)) %>%
    distinct(gene_symbol, .keep_all = TRUE)
  
  # Save all results to list
  results_list[[tp]] <- resLFC_df_anno
  
  # Write full results to file
  write_csv(resLFC_df_anno, paste0('DESeq2_tables/DESeq2_results_', tp, '_all_genes.csv'))
  
  # Filter significant genes (padj <= 0.05 and abs(log2FoldChange) > 1)
  res_sig <- resLFC_df_anno %>%
    filter(padj <= 0.05, abs(log2FoldChange) > 1)
  res_up <- res_sig %>%
    filter(log2FoldChange > 1)
  res_down <- res_sig %>%
    filter(log2FoldChange < -1)
  
  # Write significant results to files
  write_csv(res_sig, paste0('DESeq2_tables/DESeq2_results_', tp, '_significant_genes.csv'))
  write_csv(res_up, paste0('DESeq2_tables/DESeq2_results_', tp, '_upregulated_genes.csv'))
  write_csv(res_down, paste0('DESeq2_tables/DESeq2_results_', tp, '_downregulated_genes.csv'))
}

#### CDK11i/TOPKi DESeq2
setwd('../../../repos/PRMT-APA/Scripts/scripts_for_figures_in_paper/')

# Read in sample metadata
sample_info <- read_csv('CDK11i_deseq_sample_info.csv') %>%
  mutate(Treatment = factor(Treatment))

# Create vector of Salmon quant.sf file paths and assign sample names
sample_files <- paste0(sample_info$Treatment, sample_info$Replicate, '/quant.sf')
names(sample_files) <- paste0(sample_info$Treatment, '_', sample_info$Replicate)

# Import transcript-level estimates and summarize to gene-level counts
txi <- tximport(files = sample_files,
                type = 'salmon',
                tx2gene = tx2gene_map,
                countsFromAbundance = "no",
                ignoreTxVersion = TRUE)

# Round counts to integers (DESeq2 requires integer counts)
txi$counts <- round(txi$counts)

# Create DESeqDataSet object
dds_cdk11i <- DESeqDataSetFromMatrix(countData = txi$counts,
                                     colData = sample_info,
                                     design = ~ Treatment)

# Relevel Treatment factor to set control as reference
dds_cdk11i$Treatment <- relevel(dds_cdk11i$Treatment, ref = "DMSO_")

# Identify samples for each condition
condition_levels <- levels(dds_cdk11i$Treatment)
condition_A <- condition_levels[1]  # Reference level (e.g., "DMSO")
condition_B <- condition_levels[2]  # Treatment level (e.g., "CDK11i")

samples_condition_A <- colData(dds_cdk11i)$Treatment == condition_A
samples_condition_B <- colData(dds_cdk11i)$Treatment == condition_B

# Perform filtering
keep_cdk11i <- (rowSums(counts(dds_cdk11i)[, samples_condition_A] >= 5) >= 2) |
  (rowSums(counts(dds_cdk11i)[, samples_condition_B] >= 5) >= 2)
dds_cdk11i <- dds_cdk11i[keep_cdk11i, ]

# Run DESeq normalization and estimation
dds_cdk11i <- DESeq(dds_cdk11i)

# Perform variance stabilizing transformation for PCA
vsd_cdk11i <- vst(dds_cdk11i, blind = FALSE)

# PCA plot
pca_data_cdk11i <- plotPCA(vsd_cdk11i, intgroup = 'Treatment', returnData = TRUE)
percentVar_cdk11i <- round(100 * attr(pca_data_cdk11i, "percentVar"))

# Plot PCA using ggplot2
ggplot(pca_data_cdk11i, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar_cdk11i[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_cdk11i[2], "% variance")) +
  theme(text = element_text(size = 16))

# Differential expression analysis
# Get results with adjusted p-values
res_cdk11i <- results(dds_cdk11i, contrast = c("Treatment", "CDK11i_", "DMSO_"), alpha = 0.05)

# Shrink LFC estimates using apeglm
resLFC_cdk11i <- lfcShrink(dds_cdk11i, coef = "Treatment_CDK11i__vs_DMSO_", type = "apeglm")

# Convert results to data frame and add gene IDs
resLFC_df_cdk11i <- as.data.frame(resLFC_cdk11i) %>%
  rownames_to_column(var = "ensgene")

# Remove version numbers from Ensembl IDs
resLFC_df_cdk11i$ensgene <- sub("\\..*", "", resLFC_df_cdk11i$ensgene)

# Annotate genes with gene symbols
gene_annotations_cdk11i <- select(org.Hs.eg.db, keys = resLFC_df_cdk11i$ensgene, columns = "SYMBOL", keytype = "ENSEMBL")
colnames(gene_annotations_cdk11i) <- c("ensgene", "gene_symbol")
resLFC_df_anno_cdk11i <- left_join(resLFC_df_cdk11i, gene_annotations_cdk11i, by = "ensgene")

# Keep only unique genes
resLFC_df_anno_cdk11i <- resLFC_df_anno_cdk11i %>%
  group_by(ensgene) %>%
  dplyr::slice(1) %>%
  ungroup()

# Save all results
write_csv(resLFC_df_anno_cdk11i, 'DESeq2_tables/DESeq2_results_CDK11i_all_genes.csv')

# Filter significant genes (padj <= 0.05 and abs(log2FoldChange) > 1)
res_sig_cdk11i <- resLFC_df_anno_cdk11i %>%
  filter(padj <= 0.05, abs(log2FoldChange) > 1)
res_up_cdk11i <- res_sig_cdk11i %>%
  filter(log2FoldChange > 1)
res_down_cdk11i <- res_sig_cdk11i %>%
  filter(log2FoldChange < -1)

# Write significant results to files
write_csv(res_sig_cdk11i, 'DESeq2_tables/DESeq2_results_CDK11i_significant_genes.csv')
write_csv(res_up_cdk11i, 'DESeq2_tables/DESeq2_results_CDK11i_upregulated_genes.csv')
write_csv(res_down_cdk11i, 'DESeq2_tables/DESeq2_results_CDK11i_downregulated_genes.csv')

# Store results in a list for later use
results_list_cdk11i <- list(
  resLFC_df_anno = resLFC_df_anno_cdk11i,
  res_sig = res_sig_cdk11i
)

#### heatmap

# Create lists to store significant genes and all genes for each condition
sig_genes_list <- list()
all_genes_list <- list()

# Collect results from the timecourse DESeq2 analysis
for (tp in time_points) {
  # Significant genes for each time point
  sig_genes_list[[paste0("sig_", tp, "_genes")]] <- results_list[[tp]] %>%
    dplyr::filter(padj <= 0.05, abs(log2FoldChange) > 1) %>%
    dplyr::select(gene_symbol, log2FoldChange)
  
  # All genes for each time point
  all_genes_list[[paste0("all_", tp, "_genes")]] <- results_list[[tp]] %>%
    dplyr::select(gene_symbol, log2FoldChange)
}

# Collect results from the CDK11i DESeq2 analysis
sig_genes_list[["sig_CDK11i_genes"]] <- results_list_cdk11i$res_sig %>%
  dplyr::select(gene_symbol, log2FoldChange)

all_genes_list[["all_CDK11i_genes"]] <- results_list_cdk11i$resLFC_df_anno %>%
  dplyr::select(gene_symbol, log2FoldChange)

# Function to rename the 'log2FoldChange' column to include the time point
rename_lfc_column <- function(df, name_prefix) {
  colnames(df)[colnames(df) == "log2FoldChange"] <- paste0("log2FoldChange_", name_prefix)
  return(df)
}

# Rename columns to include time points or conditions
sig_genes_list <- mapply(rename_lfc_column, sig_genes_list, names(sig_genes_list), SIMPLIFY = FALSE)
all_genes_list <- mapply(rename_lfc_column, all_genes_list, names(all_genes_list), SIMPLIFY = FALSE)

# Combine significant genes from all conditions
sig_gene_symbols <- unique(unlist(lapply(sig_genes_list, function(df) df$gene_symbol)))

# Remove any rows with NA gene_symbol before merging and only keep one row per gene
all_genes_list <- lapply(all_genes_list, function(df) {
  df = df %>% filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>% dplyr::slice(1) %>%
    ungroup()
})

# Merge all gene data frames on 'gene_symbol'
all_gene_data <- purrr::reduce(
  all_genes_list,
  full_join,
  by = "gene_symbol"
)

# Collect significant gene symbols from all conditions
sig_gene_symbols <- unique(unlist(lapply(sig_genes_list, function(df) df$gene_symbol)))

# Filter to include only genes significant in any condition
sig_gene_data <- all_gene_data %>%
  filter(gene_symbol %in% sig_gene_symbols) %>%
  na.omit()  # Remove rows with NA values

# Convert tibble to data frame
sig_gene_data <- as.data.frame(sig_gene_data)

# Set 'gene_symbol' as rownames and remove it from columns
rownames(sig_gene_data) <- sig_gene_data$gene_symbol
sig_gene_data$gene_symbol <- NULL  # Remove the gene_symbol column

# Adjust column order if necessary
sig_gene_data <- sig_gene_data[, c("log2FoldChange_all_24hrs_genes",
                                   "log2FoldChange_all_48hrs_genes",
                                   "log2FoldChange_all_72hrs_genes",
                                   "log2FoldChange_all_96hrs_genes",
                                   "log2FoldChange_all_CDK11i_genes")]

# Optional: Rename columns for clarity
colnames(sig_gene_data) <- c("DMAi 24hrs", "DMAi 48hrs", "DMAi 72hrs", "DMAi 96hrs", "CDK11i/TOPKi 24hrs")

# Generate the heatmap
# Define color palette and breaks
RdBu <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(30)
logchange_breaks <- seq(-2.5, 2.5, length.out = 31)

# Create the heatmap
pheatmap(sig_gene_data,
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         color = RdBu,
         breaks = logchange_breaks,
         fontsize_row = 0.5,
         fontsize_col = 12,
         angle_col = 45,
         cutree_rows = 3)

# Extract clusters from the heatmap
# Compute the distance matrix
dist_rows <- dist(sig_gene_data, method = "manhattan")

# Perform hierarchical clustering
hc_rows <- hclust(dist_rows, method = "complete")

# Cut the tree into 3 clusters
clusters <- cutree(hc_rows, k = 3)

# Create a data frame with gene symbols and cluster assignments
gene_clusters <- data.frame(Gene = rownames(sig_gene_data), Cluster = clusters)

# Add annotation to the heatmap
annotation_row <- data.frame(Cluster = factor(clusters))
rownames(annotation_row) <- rownames(sig_gene_data)

# Generate the heatmap with annotation
pheatmap(sig_gene_data,
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         color = RdBu,
         breaks = logchange_breaks,
         fontsize_row = 0.4,
         fontsize_col = 12,
         angle_col = 45,
         cutree_rows = 3,
         annotation_row = annotation_row)


#setwd and write to file
cluster_3_genes = gene_clusters %>% filter(Cluster == 3)

write_csv(cluster_3_genes,'cluster_3_genes.csv')
write_csv(all_gene_data,'background_genes.csv')


