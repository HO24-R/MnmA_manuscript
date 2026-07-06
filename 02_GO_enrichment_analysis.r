# ============================================================
# GO ENRICHMENT DOTPLOT USING ClusterProfiler + org.EcK12.eg.db
#
# Purpose:
# Performs separate GO enrichment analyses for upregulated and downregulated genes

#
# Note: The list of genes here have been curated based on log2FC of <= -1.50 and >= 1.50; padj = 0.05
# No additional differential-expression cutoffs are applied in this script. All genes present in the input files are used
# after validation against org.EcK12.eg.db.

# ============================================================
# 1. CLEAN ENVIRONMENT AND LOAD REQUIRED PACKAGES
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(clusterProfiler)
  library(org.EcK12.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})


# ============================================================
# 2. USER SETTINGS
#
# Define input files, output directory, gene identifier column,
# GO ontology, number of terms displayed, and output format.
# ============================================================

UP_FILE   <- "Upregulated_***.csv"
DOWN_FILE <- "Downregulated_***.csv"

OUTDIR <- "Name_of_folder_to_create"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

GENE_COL <- "GeneName"       # Column containing gene symbols
ONT <- "BP"                  # BP / MF / CC
SHOW_N_PER_DIR <- 12         # Top terms shown per direction
SIMPLIFY_TERMS <- TRUE       # Remove redundant GO terms
SAVE_PDF <- TRUE             # Save PDF output in addition to PNG


# ============================================================
# 3. HELPER FUNCTIONS
#
# Functions for:
#   - converting GeneRatio to numeric format;
#   - reading and validating gene lists;
#   - running GO enrichment analysis.
# ============================================================

# Convert GeneRatio from character format (e.g., "10/100") to a numeric proportion (e.g., 0.10).

parse_generatio <- function(gr) {
  
  parts <- str_split(gr, "/", simplify = TRUE)
  
  as.numeric(parts[, 1]) /
    as.numeric(parts[, 2])
}


# Read gene lists, retain valid E. coli gene symbols,
# and save mapped and unmapped gene lists.

read_gene_list <- function(path, label) {
  
  df <- read_csv(path, show_col_types = FALSE)
  
# Confirm that the required gene column is present.
  
  if (!(GENE_COL %in% names(df))) {
    
    stop(
      "Column '", GENE_COL,
      "' not found in: ", path,
      "\nColumns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
# Extract unique gene symbols.
  
  genes_raw <- df %>%
    transmute(
      gene = as.character(.data[[GENE_COL]])
    ) %>%
    filter(
      !is.na(gene),
      gene != ""
    ) %>%
    distinct() %>%
    pull(gene)
  
  # Validate gene symbols against the E. coli annotation database.
  
  valid_syms <- AnnotationDbi::keys(
    org.EcK12.eg.db,
    keytype = "SYMBOL"
  )
  
  genes <- intersect(
    genes_raw,
    valid_syms
  )
  
  not_mapped <- setdiff(
    genes_raw,
    genes
  )
  
  # Save mapped and unmapped gene lists for QC.
  
  writeLines(
    genes,
    file.path(
      OUTDIR,
      paste0(label, "_genes_used_SYMBOL.txt")
    )
  )
  
  
  if (length(not_mapped) > 0) {
    
    writeLines(
      not_mapped,
      file.path(
        OUTDIR,
        paste0(label, "_genes_not_mapped.txt")
      )
    )
  }
  
# Report gene-mapping summary.
  
  message(
    label,
    ": genes in file = ", length(genes_raw),
    " | mapped SYMBOL = ", length(genes),
    " | not mapped = ", length(not_mapped)
  )
  
  
  if (length(genes) < 5) {
    
    stop(
      "Too few mapped genes for ",
      label,
      ". Check GeneName format."
    )
  }
  
  
  genes
}


# Perform GO enrichment analysis and optionally remove semantically redundant GO terms.

run_enrichGO <- function(genes, direction_label) {
  
  ego <- enrichGO(
    gene          = genes,
    OrgDb         = org.EcK12.eg.db,
    keyType       = "SYMBOL",
    ont           = ONT,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  
  
  if (SIMPLIFY_TERMS) {
    
    ego <- suppressWarnings(
      simplify(
        ego,
        cutoff = 0.7,
        by = "p.adjust",
        select_fun = min
      )
    )
  }
  
  
  ego_df <- as.data.frame(ego)
  
  ego_df$Direction <- direction_label
  
  ego_df
}


# ============================================================
# 4. LOAD AND VALIDATE UPREGULATED AND DOWNREGULATED GENES
# ============================================================

up_genes <- read_gene_list(
  UP_FILE,
  "Upregulated"
)

down_genes <- read_gene_list(
  DOWN_FILE,
  "Downregulated"
)

# ============================================================
# 5. RUN GO ENRICHMENT SEPARATELY BY DIRECTION
#
# Upregulated and downregulated genes are analyzed separately so their enriched biological processes can be compared.

up_df <- run_enrichGO(
  up_genes,
  "Upregulated"
)

down_df <- run_enrichGO(
  down_genes,
  "Downregulated"
)

# Save complete GO enrichment tables.

write_csv(
  up_df,
  file.path(
    OUTDIR,
    paste0(
      "enrichGO_",
      ONT,
      "_Upregulated_NO_CUTOFF.csv"
    )
  )
)

write_csv(
  down_df,
  file.path(
    OUTDIR,
    paste0(
      "enrichGO_",
      ONT,
      "_Downregulated_NO_CUTOFF.csv"
    )
  )
)


# Stop if no GO terms are returned for either gene list.

if (nrow(up_df) == 0 && nrow(down_df) == 0) {
  
  stop(
    "No GO terms returned for either list. ",
    "Check mapping / OrgDb / ontology."
  )
}

# ============================================================
# 6. SELECT TOP GO TERMS FOR VISUALIZATION
#
# Terms are ranked by adjusted p-value, and the top terms from each direction are retained for plotting.
# ============================================================

top_up <- up_df %>%
  filter(!is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = SHOW_N_PER_DIR)


top_down <- down_df %>%
  filter(!is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = SHOW_N_PER_DIR)


# Combine both directions and prepare GeneRatio and
# GO-term descriptions for plotting.

plot_df <- bind_rows(
  top_down,
  top_up
) %>%
  mutate(
    GeneRatio_num = parse_generatio(GeneRatio),
    
    Description_wrapped = str_wrap(
      Description,
      width = 45
    )
  )

# ============================================================
# 7. ORDER GO TERMS ON THE SHARED Y-AXIS
#
# GO terms are ordered using the best adjusted p-value observed across the two directions, with stronger terms placed higher.
# ============================================================

term_order <- plot_df %>%
  group_by(Description_wrapped) %>%
  summarise(
    best_p = min(
      p.adjust,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(best_p) %>%
  pull(Description_wrapped)


plot_df <- plot_df %>%
  mutate(
    
    Direction = factor(
      Direction,
      levels = c(
        "Downregulated",
        "Upregulated"
      )
    ),
    
    Description_wrapped = factor(
      Description_wrapped,
      levels = rev(term_order)
    )
  )

# ============================================================
# 8. GENERATE TWO-COLUMN GO DOTPLOT
#
# Dot size represents GeneRatio, while dot color represents
# the adjusted p-value of each enriched GO term.
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = Direction,
    y = Description_wrapped
  )
) +
  
  geom_point(
    aes(
      size = GeneRatio_num,
      color = p.adjust
    )
  ) +
  
  scale_size_continuous(
    name = "GeneRatio"
  ) +
  
  scale_color_gradient(
    name = "p.adjust",
    low = "blue",
    high = "red",
    trans = "reverse"
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid.minor = element_blank(),
    
    axis.text.y = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    )
  )


print(p)


# ============================================================
# 9. SAVE GO DOTPLOT
#
# The figure is saved as a high-resolution PNG and, optionally, as a publication-quality PDF.
# ============================================================

ggsave(
  file.path(
    OUTDIR,
    paste0(
      "GO_dotplot_",
      ONT,
      "_Up_vs_Down.png"
    )
  ),
  p,
  width = 10,
  height = 10,
  dpi = 300
)


if (SAVE_PDF) {
  
  ggsave(
    file.path(
      OUTDIR,
      paste0(
        "GO_dotplot_",
        ONT,
        "_Up_vs_Down.pdf"
      )
    ),
    p,
    width = 10,
    height = 10
  )
}


# ============================================================
# 10. COMPLETION MESSAGE
# ============================================================

message(
  "Done. Outputs written to: ",
  OUTDIR
)

# ============================================================
# END OF SCRIPT; ENJOY
# ============================================================




#################################################### TO USE THE SINGLE DESEQ OUTPUT FILE FOR GO ANALYSIS ################

# ============================================================
# GO ENRICHMENT ANALYSIS: ONE DESEQ2 FILE
#
# Purpose:
# Reads one DESeq2 differential-expression results file,
# separates upregulated and downregulated genes using a
# log2 fold-change cutoff, performs GO enrichment analysis
# separately for each direction, and generates a two-column
# GO dotplot.
#
# Method:
# clusterProfiler + org.EcK12.eg.db
#
# Filtering:
# Genes are selected using log2 fold-change only:
#   Upregulated:   log2FC >= 1.5
#   Downregulated: log2FC <= -1.5
#
# Outputs:
#   - GO enrichment tables for each direction
#   - Mapped and unmapped gene lists
#   - Plotting data table
#   - GO dotplot in PDF, SVG, PNG, and TIFF formats
# ============================================================


# ============================================================
# 1. USER SETTINGS
#
# Define the DESeq2 input file, output directory, gene and
# log2FC columns, GO ontology, and plotting parameters.
# ============================================================

DE_FILE <- "file_name.csv"

OUTDIR <- "name_of_folder"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

GENE_COL <- "GeneName"
LFC_COL  <- "log2FoldChange"

ONT <- "BP"              # BP / MF / CC
LFC_CUTOFF <- 1.5        # Filter genes using log2FC only
SHOW_N_PER_DIR <- 12     # Top GO terms shown per direction
SIMPLIFY_TERMS <- TRUE   # Remove redundant GO terms


# ============================================================
# 2. HELPER FUNCTIONS
#
# Functions for:
#   - converting GeneRatio to numeric format;
#   - obtaining valid E. coli gene symbols;
#   - performing GO enrichment analysis.
# ============================================================


# Convert GeneRatio from character format, e.g. "10/100",
# to a numeric proportion, e.g. 0.10.

parse_generatio <- function(gr) {
  parts <- stringr::str_split(gr, "/", simplify = TRUE)
  as.numeric(parts[, 1]) / as.numeric(parts[, 2])
}


# Retrieve valid E. coli gene symbols from the annotation database.

get_valid_symbols <- function() {
  AnnotationDbi::keys(org.EcK12.eg.db, keytype = "SYMBOL")
}


# Perform GO enrichment analysis and optionally remove
# semantically redundant GO terms.

run_enrichGO <- function(genes, direction_label) {
  
  ego <- clusterProfiler::enrichGO(
    gene          = unique(genes),
    OrgDb         = org.EcK12.eg.db,
    keyType       = "SYMBOL",
    ont           = ONT,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(data.frame())
  }
  
  if (SIMPLIFY_TERMS) {
    ego <- suppressWarnings(
      clusterProfiler::simplify(
        ego,
        cutoff = 0.7,
        by = "p.adjust",
        select_fun = min
      )
    )
  }
  
  ego_df <- as.data.frame(ego)
  ego_df$Direction <- direction_label
  
  ego_df
}


# ============================================================
# 3. LOAD AND PREPARE DESEQ2 RESULTS
#
# Read the DESeq2 results file, verify the required columns,
# and retain genes with valid gene names and log2FC values.
# ============================================================

de_df <- readr::read_csv(DE_FILE, show_col_types = FALSE)

required_cols <- c(GENE_COL, LFC_COL)

missing_cols <- setdiff(required_cols, names(de_df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    "\nColumns found: ",
    paste(names(de_df), collapse = ", ")
  )
}

de_df <- de_df %>%
  dplyr::mutate(
    GeneName = as.character(.data[[GENE_COL]]),
    log2FoldChange = as.numeric(.data[[LFC_COL]])
  ) %>%
  dplyr::filter(
    !is.na(GeneName),
    GeneName != "",
    !is.na(log2FoldChange)
  )


# ============================================================
# 4. IDENTIFY UPREGULATED AND DOWNREGULATED GENES
#
# Split genes by direction using only the specified log2FC
# cutoff and retain symbols mapped in org.EcK12.eg.db.
# ============================================================

valid_syms <- get_valid_symbols()

up_genes_raw <- de_df %>%
  dplyr::filter(log2FoldChange >= LFC_CUTOFF) %>%
  dplyr::pull(GeneName) %>%
  unique()

down_genes_raw <- de_df %>%
  dplyr::filter(log2FoldChange <= -LFC_CUTOFF) %>%
  dplyr::pull(GeneName) %>%
  unique()

up_genes <- intersect(up_genes_raw, valid_syms)
down_genes <- intersect(down_genes_raw, valid_syms)

up_not_mapped <- setdiff(up_genes_raw, up_genes)
down_not_mapped <- setdiff(down_genes_raw, down_genes)


# Save mapped gene lists used for GO enrichment.

writeLines(
  up_genes,
  file.path(OUTDIR, "Upregulated_genes_used_SYMBOL.txt")
)

writeLines(
  down_genes,
  file.path(OUTDIR, "Downregulated_genes_used_SYMBOL.txt")
)


# Save genes that could not be mapped.

if (length(up_not_mapped) > 0) {
  writeLines(
    up_not_mapped,
    file.path(OUTDIR, "Upregulated_genes_not_mapped.txt")
  )
}

if (length(down_not_mapped) > 0) {
  writeLines(
    down_not_mapped,
    file.path(OUTDIR, "Downregulated_genes_not_mapped.txt")
  )
}


# Report gene-selection and mapping summaries.

message(
  "Upregulated: genes passing log2FC cutoff = ", length(up_genes_raw),
  " | mapped SYMBOL = ", length(up_genes),
  " | not mapped = ", length(up_not_mapped)
)

message(
  "Downregulated: genes passing log2FC cutoff = ", length(down_genes_raw),
  " | mapped SYMBOL = ", length(down_genes),
  " | not mapped = ", length(down_not_mapped)
)


# Warn if few genes are available for enrichment analysis.

if (length(up_genes) < 5) {
  warning(
    "Fewer than 5 mapped upregulated genes. ",
    "Upregulated GO may be weak or empty."
  )
}

if (length(down_genes) < 5) {
  warning(
    "Fewer than 5 mapped downregulated genes. ",
    "Downregulated GO may be weak or empty."
  )
}


# ============================================================
# 5. RUN GO ENRICHMENT SEPARATELY BY DIRECTION
#
# Upregulated and downregulated genes are analyzed separately
# to identify direction-specific enriched biological processes.
# ============================================================

up_df <- run_enrichGO(
  up_genes,
  "Upregulated"
)

down_df <- run_enrichGO(
  down_genes,
  "Downregulated"
)


# Save complete GO enrichment results.

readr::write_csv(
  up_df,
  file.path(
    OUTDIR,
    paste0(
      "enrichGO_",
      ONT,
      "_Upregulated_LFC_",
      LFC_CUTOFF,
      ".csv"
    )
  )
)

readr::write_csv(
  down_df,
  file.path(
    OUTDIR,
    paste0(
      "enrichGO_",
      ONT,
      "_Downregulated_LFC_",
      LFC_CUTOFF,
      ".csv"
    )
  )
)


if (nrow(up_df) == 0 && nrow(down_df) == 0) {
  stop(
    "No GO terms returned for either Upregulated ",
    "or Downregulated genes."
  )
}


# ============================================================
# 6. SELECT TOP GO TERMS FOR VISUALIZATION
#
# Terms are ranked by adjusted p-value, and the top terms from
# each direction are retained for plotting.
# ============================================================

top_up <- up_df %>%
  dplyr::filter(!is.na(p.adjust)) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = SHOW_N_PER_DIR)

top_down <- down_df %>%
  dplyr::filter(!is.na(p.adjust)) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = SHOW_N_PER_DIR)


# Combine selected terms and prepare GeneRatio and GO-term
# descriptions for plotting.

plot_df <- dplyr::bind_rows(
  top_down,
  top_up
) %>%
  dplyr::mutate(
    GeneRatio_num = parse_generatio(GeneRatio),
    Description_wrapped = stringr::str_wrap(
      Description,
      width = 45
    )
  )


if (nrow(plot_df) == 0) {
  stop("No GO terms available for plotting.")
}


# ============================================================
# 7. ORDER GO TERMS ON THE SHARED Y-AXIS
#
# GO terms are ordered using the best adjusted p-value observed
# across both directions, with stronger terms positioned higher.
# ============================================================

term_order <- plot_df %>%
  dplyr::group_by(Description_wrapped) %>%
  dplyr::summarise(
    best_p = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(best_p) %>%
  dplyr::pull(Description_wrapped)


plot_df <- plot_df %>%
  dplyr::mutate(
    Direction = factor(
      Direction,
      levels = c(
        "Downregulated",
        "Upregulated"
      )
    ),
    Description_wrapped = factor(
      Description_wrapped,
      levels = rev(term_order)
    )
  )


# Save the exact data used to generate the dotplot.

readr::write_csv(
  plot_df,
  file.path(
    OUTDIR,
    paste0(
      "GO_dotplot_",
      ONT,
      "_Up_vs_Down_plot_data.csv"
    )
  )
)


# ============================================================
# 8. GENERATE TWO-COLUMN GO DOTPLOT
#
# Dot size represents GeneRatio, while dot color represents
# the adjusted p-value of each enriched GO term.
# ============================================================

p <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = Direction,
    y = Description_wrapped
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = GeneRatio_num,
      color = p.adjust
    ),
    alpha = 0.95
  ) +
  ggplot2::scale_size_continuous(
    name = "GeneRatio"
  ) +
  ggplot2::scale_color_gradient(
    name = "p.adjust",
    low = "blue",
    high = "red",
    trans = "reverse"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(
    base_size = 11
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    legend.text = ggplot2::element_text(
      size = 10
    )
  )

print(p)


# ============================================================
# 9. SAVE GO DOTPLOT
#
# Save the final figure in vector and high-resolution raster
# formats for manuscript preparation and downstream use.
# ============================================================

PLOT_NAME <- file.path(
  OUTDIR,
  paste0(
    "GO_dotplot_",
    ONT,
    "_Up_vs_Down_LFC_",
    LFC_CUTOFF
  )
)


# PDF

ggplot2::ggsave(
  paste0(PLOT_NAME, ".pdf"),
  p,
  width = 8.2,
  height = 5.0
)


# SVG

ggplot2::ggsave(
  paste0(PLOT_NAME, ".svg"),
  p,
  width = 8.2,
  height = 5.0,
  device = svglite::svglite
)


# PNG

ggplot2::ggsave(
  paste0(PLOT_NAME, ".png"),
  p,
  width = 8.2,
  height = 5.0,
  dpi = 600
)


# TIFF

ggplot2::ggsave(
  paste0(PLOT_NAME, ".tiff"),
  p,
  width = 8.2,
  height = 5.0,
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# 10. COMPLETION MESSAGE
# ============================================================

message(
  "Done. Outputs written to: ",
  OUTDIR
)


# ============================================================
# END OF SCRIPT
# ============================================================