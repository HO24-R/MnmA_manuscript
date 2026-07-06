# ============================================================
# GO ENRICHMENT DOTPLOT USING ClusterProfiler + org.EcK12.eg.db
#
# Purpose:
# Performs separate GO enrichment analyses for upregulated and downregulated genes

#
# Note:
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
# END OF SCRIPT
# ============================================================