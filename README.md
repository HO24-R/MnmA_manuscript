# MnmA Manuscript Code

This repository contains R scripts used for RNA-seq differential expression analysis and downstream gene annotation for the MnmA tRNA modification manuscript.

## Script included

### 01_DESeq2_WT_treated_vs_WT_untreated.R

This script performs DESeq2 differential expression analysis comparing WT treated versus WT untreated samples.

## Required inputs

Place STAR gene count files here:

```text
Data/Raw/ (each file should end with "ReadsPerGene.out.tab)

The genome annotation file is placed in: Data/Annotation/

Results are saved in: Results/ (main output files include "WT_treated_vs_WT_untreated_DE.csv"
and "WT_treated_vs_WT_untreated_DE_with_Gene_Names.csv"

Required packages include: 
DESeq2
dplyr
rtracklayer
