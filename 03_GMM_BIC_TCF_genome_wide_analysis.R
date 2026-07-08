# ============================================================
# GENOME-WIDE CODON FREQUENCY CLUSTERING AND CODON BIAS ANALYSIS
#
# Purpose:
#   1. Analyze gene-level codon frequency data.
#   2. Cluster genes using Gaussian Mixture Models (GMM).
#   3. Compare A/U-ending and G/C-ending codon preferences across clusters using bootstrap confidence intervals.
#   4. Generate an A/U vs G/C dumbbell plot.
#
# Input:
#   GCF_000184185_codon_frequencies.xlsx
#   Sheet: TotFreqs
# ============================================================


# ============================================================
# 1. CLEAN ENVIRONMENT AND LOAD PACKAGES
# ============================================================

rm(list = ls())

# CRAN packages
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
library(ggplot2)
library(mclust)
library(ggrepel)

# Bioconductor packages
library(clusterProfiler)
library(org.EcK12.eg.db)
library(enrichplot)

# Set seed for reproducibility
set.seed(333)


# ============================================================
# 2. LOAD AND NORMALIZE CODON FREQUENCY DATA
#
# Remove stop codons and methionine, Z-score each codon
# frequency across genes, and convert codon labels from
# DNA notation (T) to RNA notation (U).
# ============================================================

TCF <- read_excel("GCF_000184185_codon_frequencies.xlsx", sheet = "TotFreqs")

# Remove stop codons and methionine
TCF <- TCF %>%
  dplyr::select(-dplyr::any_of(c("STOPTAA", "STOPTAG", "STOPTGA", "MetATG")))

# Rename first column as gene identifier
names(TCF)[1] <- "Gene"

# Identify raw codon-frequency columns
codon_cols_raw <- setdiff(names(TCF), "Gene")

# Z-score each codon across genes and replace missing values with zero
TCF <- TCF %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(codon_cols_raw), ~ as.numeric(scale(.)))
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(codon_cols_raw), ~ tidyr::replace_na(., 0))
  )

# Convert T to U in the codon portion of column names
colnames(ICF)[-1] <- colnames(ICF)[-1] %>%
  str_replace("(...)([A-Z]{3})$", function(x) {
    aa <- str_match(x, "(...)([A-Z]{3})")[, 2]
    codon <- str_match(x, "(...)([A-Z]{3})")[, 3]
    paste0(aa, str_replace_all(codon, "T", "U"))
  })

# Identify codon columns after conversion to RNA notation
codon_cols <- names(TCF) %>%
  keep(~ str_detect(.x, "^[A-Za-z]{3}[AUCG]{3}$"))


# ============================================================
# 3. CLUSTER GENES USING GAUSSIAN MIXTURE MODELS
#
# Compare GMM solutions across 1–25 clusters using BIC and
# fit the best-supported model.
# ============================================================

X <- TCF %>%
  dplyr::select(dplyr::all_of(codon_cols))

# Calculate BIC across candidate cluster numbers
bic <- mclustBIC(X, G = 1:25)

# Plot BIC model comparison
plot(bic, legendArgs = list(x = "topright", cex = 0.8))

# Display best model and cluster solution
best <- summary(bic)
best

# Fit the BIC-selected GMM
fit <- Mclust(X, x = bic)
summary(fit)

# Add cluster assignments to codon-frequency table
TCF <- TCF %>%
  dplyr::mutate(
    ClusterNum = as.integer(fit$classification),
    Cluster = factor(paste0("Cluster ", ClusterNum))
  )


# ============================================================
# 4. CLASSIFY CODONS BY THIRD-POSITION NUCLEOTIDE
#
# Group codons as A/U-ending or G/C-ending according to the
# nucleotide at the third codon position.
# ============================================================

codon_ending_map <- tibble(
  Codon = codon_cols,
  Ending = dplyr::case_when(
    str_sub(Codon, -1) %in% c("A", "U") ~ "A/U-ending",
    str_sub(Codon, -1) %in% c("G", "C") ~ "G/C-ending",
    TRUE ~ NA_character_
  )
)


# ============================================================
# 5. CALCULATE GENE-LEVEL A/U AND G/C CODON SCORES
#
# Calculate the mean codon Z-score for A/U-ending and
# G/C-ending codons within each gene.
# ============================================================

gene_end_means <- TCF %>%
  dplyr::select(Gene, Cluster, dplyr::all_of(codon_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(codon_cols),
    names_to = "Codon",
    values_to = "Z"
  ) %>%
  dplyr::left_join(codon_ending_map, by = "Codon") %>%
  dplyr::filter(!is.na(Ending)) %>%
  dplyr::group_by(Gene, Cluster, Ending) %>%
  dplyr::summarise(
    GeneMeanZ = mean(Z, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 6. CALCULATE BOOTSTRAP CONFIDENCE INTERVALS
#
# Estimate the cluster-level mean Z-score and 95% bootstrap
# confidence interval for each codon-ending group.
# ============================================================

boot_mean_ci <- function(x, nboot = 2000, conf = 0.95) {
  
  x <- x[is.finite(x)]
  
  if (length(x) < 2) {
    return(tibble(
      mean = mean(x),
      lo = NA_real_,
      hi = NA_real_,
      n = length(x)
    ))
  }
  
  boots <- replicate(nboot, mean(sample(x, replace = TRUE)))
  
  tibble(
    mean = mean(x),
    lo = unname(quantile(boots, (1 - conf) / 2)),
    hi = unname(quantile(boots, 1 - (1 - conf) / 2)),
    n = length(x)
  )
}


summary_ci <- gene_end_means %>%
  dplyr::group_by(Cluster, Ending) %>%
  dplyr::summarise(
    boot_mean_ci(GeneMeanZ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Ending = factor(Ending, levels = c("A/U-ending", "G/C-ending"))
  )


# ============================================================
# 7. ORDER CLUSTERS BY A/U-ENDING CODON SCORE
#
# Arrange clusters from highest to lowest mean A/U-ending
# codon Z-score for visualization.
# ============================================================

cluster_order <- summary_ci %>%
  dplyr::filter(Ending == "A/U-ending") %>%
  dplyr::arrange(desc(mean)) %>%
  dplyr::pull(Cluster)

summary_ci <- summary_ci %>%
  dplyr::mutate(
    Cluster = factor(Cluster, levels = cluster_order)
  )


# ============================================================
# 8. GENERATE A/U VS G/C DUMBBELL PLOT
#
# Plot cluster-level mean codon Z-scores and 95% bootstrap
# confidence intervals for A/U-ending and G/C-ending codons.
# ============================================================

p_dumbbell <- ggplot(
  summary_ci,
  aes(y = Cluster, color = Ending)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 1.0
  ) +
  geom_errorbarh(
    aes(xmin = lo, xmax = hi),
    height = 0,
    linewidth = 1.2
  ) +
  geom_point(aes(x = lo), size = 1.8) +
  geom_point(aes(x = hi), size = 1.8) +
  geom_point(aes(x = mean), size = 4.0) +
  scale_color_manual(
    values = c("A/U-ending" = "gold", "G/C-ending" = "#7B3294")
  ) +
  labs(
    x = "Average Z-score (± 95% CI)",
    y = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold")
  )

print(p_dumbbell)


# ============================================================
# END OF SCRIPT
# ============================================================