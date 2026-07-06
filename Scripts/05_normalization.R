# ─────────────────────────────────────────────────────────────────────────────
# 05_normalization.R — Count Normalization
# RNA-Seq Pipeline | Karim Ghazouani
# ─────────────────────────────────────────────────────────────────────────────
# DESCRIPTION:
#   Normalizes raw counts using DESeq2 size factors (median-of-ratios method).
#   Produces normalized count tables and QC plots (boxplot, PCA) to verify
#   that normalization corrects for library size differences.
#
# INPUT:
#   results/counts/counts_matrix.csv — raw count matrix
#   data/samples.txt                 — sample metadata
#
# OUTPUT:
#   results/normalization/normalized_counts.csv — DESeq2 normalized counts
#   results/normalization/vst_counts.csv        — variance-stabilized counts
#   results/normalization/plots/                — QC plots
# ─────────────────────────────────────────────────────────────────────────────

library(DESeq2)
library(ggplot2)
library(reshape2)

# ── Parameters ────────────────────────────────────────────────────────────────
counts_file   <- "results/counts/counts_matrix.csv"
metadata_file <- "data/samples.txt"
output_dir    <- "results/normalization"
plots_dir     <- file.path(output_dir, "plots")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
counts_matrix      <- read.csv(counts_file, row.names = 1)
counts_matrix      <- as.matrix(counts_matrix)

metadata           <- read.table(metadata_file, header = TRUE, row.names = "sample")
metadata$condition <- factor(metadata$condition)

# Ensure column order matches metadata
counts_matrix <- counts_matrix[, rownames(metadata)]

# ── Build DESeqDataSet ────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData   = metadata,
                              design    = ~ condition)

# Filter lowly expressed genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes retained after filtering:", nrow(dds), "\n")

# ── Estimate size factors and normalize ───────────────────────────────────────
dds                <- estimateSizeFactors(dds)
normalized_counts  <- counts(dds, normalized = TRUE)

cat("\nSize factors:\n")
print(sizeFactors(dds))

# ── Variance-stabilizing transformation (for visualization) ──────────────────
vsd         <- vst(dds, blind = TRUE)
vst_matrix  <- assay(vsd)

# ── Save outputs ──────────────────────────────────────────────────────────────
write.csv(normalized_counts,
          file.path(output_dir, "normalized_counts.csv"))
write.csv(vst_matrix,
          file.path(output_dir, "vst_counts.csv"))

# ── QC Plot 1: Boxplot of log2 counts before and after normalization ──────────
raw_log   <- log2(counts_matrix + 1)
norm_log  <- log2(normalized_counts + 1)

plot_boxplot <- function(mat, title, filepath) {
  df <- melt(mat, varnames = c("gene", "sample"))
  p  <- ggplot(df, aes(x = sample, y = value, fill = sample)) +
    geom_boxplot(outlier.size = 0.3, show.legend = FALSE) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = NULL, y = "log2(counts + 1)")
  ggsave(filepath, p, width = 10, height = 5)
}

plot_boxplot(raw_log,
             "Raw counts (log2)",
             file.path(plots_dir, "boxplot_raw.png"))

plot_boxplot(norm_log,
             "Normalized counts (log2)",
             file.path(plots_dir, "boxplot_normalized.png"))

# ── QC Plot 2: PCA on VST data ────────────────────────────────────────────────
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  theme_minimal(base_size = 12) +
  ggtitle("PCA — Variance-Stabilized Counts")

ggsave(file.path(plots_dir, "pca_vst.png"), pca_plot, width = 7, height = 5)

cat("\nNormalization complete.\n")
cat("Outputs saved to:", output_dir, "\n")
