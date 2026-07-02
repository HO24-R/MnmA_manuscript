# Script 01: DESeq2 analysis for WT treated vs WT untreated
# Project: MnmA tRNA modification manuscript
#
# Purpose:
#   1. Build a gene count matrix from STAR ReadsPerGene.out.tab files
#   2. Run DESeq2 differential expression analysis
#   3. Map gene IDs to gene names using a GFF3 annotation file
#   4. Save DESeq2 results with gene names
#
# Inputs:
#   - data/raw/*ReadsPerGene.out.tab
#   - data/annotation/Annot.gff3
#
# Outputs:
#   - results/WT_treated_vs_WT_untreated_DE.csv
#   - results/WT_treated_vs_WT_untreated_DE_with_Gene_Names.csv
#
# Author: Humphrey Omeoga
# ============================================================


# -----------------------------
# 1. Clear environment
# -----------------------------
rm(list = ls())


# -----------------------------
# 2. Load required packages
# -----------------------------
library(DESeq2)
library(dplyr)
library(rtracklayer)


# -----------------------------
# 3. Define input and output paths
# -----------------------------
counts_dir <- "data/raw"
annotation_file <- "data/annotation/Annot.gff3"
outdir <- "results"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# -----------------------------
# 4. Load STAR ReadsPerGene files
# -----------------------------
file.list <- list.files(
  path = counts_dir,
  pattern = "ReadsPerGene.out.tab$",
  full.names = TRUE
)

if (length(file.list) == 0) {
  stop("No ReadsPerGene.out.tab files found in data/raw/")
}

counts.files <- lapply(file.list, read.table, skip = 4)


# -----------------------------
# 5. Build count matrix
# -----------------------------
counts <- as.data.frame(
  sapply(counts.files, function(x) x[, 2])
)

colnames(counts) <- basename(file.list)
rownames(counts) <- counts.files[[1]]$V1


# -----------------------------
# 6. Define sample conditions
# IMPORTANT:
# This assumes the first 4 files are WT_Treated
# and the next 4 files are WT_Untreated.
# Make sure file.list is in the correct order.
# -----------------------------
condition <- c(
  rep("WT_Treated", 4),
  rep("WT_Untreated", 4)
)

if (length(condition) != ncol(counts)) {
  stop("Number of conditions does not match number of count files.")
}

sampleTable <- data.frame(
  sampleName = colnames(counts),
  condition = condition
)

rownames(sampleTable) <- sampleTable$sampleName


# -----------------------------
# 7. Run DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sampleTable,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "WT_Treated", "WT_Untreated")
)

res_ordered <- res[order(res$padj), ]


# -----------------------------
# 8. Save DESeq2 results
# -----------------------------
res_df <- as.data.frame(res_ordered)
res_df$genes <- rownames(res_df)

write.csv(
  res_df,
  file = file.path(outdir, "WT_treated_vs_WT_untreated_DE.csv"),
  row.names = FALSE
)


# -----------------------------
# 9. Import GFF3 annotation
# -----------------------------
gff_data <- import.gff3(annotation_file, format = "GFF3")

gene_IDs <- mcols(gff_data)$gene_id
gene_names <- mcols(gff_data)$Name

gene_mapping <- data.frame(
  genes = gene_IDs,
  GeneName = gene_names,
  stringsAsFactors = FALSE
)

gene_mapping <- gene_mapping %>%
  filter(!is.na(genes)) %>%
  distinct(genes, .keep_all = TRUE)


# -----------------------------
# 10. Add gene names to DESeq2 table
# -----------------------------
merged_data <- left_join(
  res_df,
  gene_mapping,
  by = "genes"
)

write.csv(
  merged_data,
  file = file.path(outdir, "WT_treated_vs_WT_untreated_DE_with_Gene_Names.csv"),
  row.names = FALSE
)


# -----------------------------
# End of script
# -----------------------------