# ─────────────────────────────────────────────────────────────────────────────
# 06_differential_expression.R — Differential Expression Analysis
# RNA-Seq Pipeline | Dr. Hicham Charoute
# ─────────────────────────────────────────────────────────────────────────────
# DESCRIPTION:
#   Runs DESeq2 differential expression analysis (Wald test).
#   Annotates results with HGNC gene symbols.
#   Produces MA plot, volcano plot, and heatmap of top DE genes.
#
# INPUT:
#   results/counts/counts_matrix.csv — raw count matrix
#   data/samples.txt                 — sample metadata
#
# OUTPUT:
#   results/de/de_results_all.csv        — all genes with stats
#   results/de/de_results_significant.csv — significant DE genes (padj < 0.05)
#   results/de/plots/                    — MA plot, volcano, heatmap
# ─────────────────────────────────────────────────────────────────────────────

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(dplyr)

# ── Parameters ────────────────────────────────────────────────────────────────
counts_file         <- "results/counts/counts_matrix.csv"
metadata_file       <- "data/samples.txt"
output_dir          <- "results/de"
plots_dir           <- file.path(output_dir, "plots")
reference_condition <- "control"
padj_threshold      <- 0.05
lfc_threshold       <- 1.5
top_n_heatmap       <- 50

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
counts_matrix      <- as.matrix(read.csv(counts_file, row.names = 1))
metadata           <- read.table(metadata_file, header = TRUE, row.names = "sample")
metadata$condition <- factor(metadata$condition)
counts_matrix      <- counts_matrix[, rownames(metadata)]

# ── DESeq2 ────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData   = metadata,
                              design    = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
dds$condition <- relevel(dds$condition, ref = reference_condition)

dds <- DESeq(dds)
res <- results(dds, alpha = padj_threshold)

cat("\n── DESeq2 Results Summary ──────────────────────────────────────────\n")
summary(res)

# ── Annotation ────────────────────────────────────────────────────────────────
ensembl_ids <- rownames(res)

symbols <- mapIds(org.Hs.eg.db,
                  keys      = ensembl_ids,
                  keytype   = "ENSEMBL",
                  column    = "SYMBOL",
                  multiVals = "first")

symbols[is.na(symbols)] <- ensembl_ids[is.na(symbols)]
res$genename            <- make.unique(symbols)

# ── Save full and significant results ─────────────────────────────────────────
res_df  <- as.data.frame(res)
res_sig <- res_df[!is.na(res_df$padj) & res_df$padj < padj_threshold, ]
res_sig <- res_sig[order(res_sig$padj), ]

write.csv(res_df,  file.path(output_dir, "de_results_all.csv"))
write.csv(res_sig, file.path(output_dir, "de_results_significant.csv"))

cat("\nSignificant DE genes (padj <", padj_threshold, "):", nrow(res_sig), "\n")

# ── Plot 1: MA Plot ───────────────────────────────────────────────────────────
png(file.path(plots_dir, "ma_plot.png"), width = 900, height = 600, res = 120)
plotMA(res, alpha = padj_threshold, main = "MA Plot", ylim = c(-5, 5))
dev.off()

# ── Plot 2: Volcano Plot ──────────────────────────────────────────────────────
volcano_data <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(regulation = case_when(
    padj < padj_threshold & log2FoldChange >  lfc_threshold ~ "Upregulated",
    padj < padj_threshold & log2FoldChange < -lfc_threshold ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), alpha = 0.6, size = 1.4) +
  scale_color_manual(values = c("Upregulated"   = "#E41A1C",
                                "Downregulated" = "#377EB8",
                                "Not significant" = "grey70")) +
  geom_text_repel(
    data         = subset(volcano_data, padj < 0.01 & abs(log2FoldChange) > 2),
    aes(label    = genename),
    max.overlaps = 20,
    size         = 3) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
             linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(padj_threshold),
             linetype = "dashed", alpha = 0.5) +
  theme_minimal(base_size = 12) +
  labs(title = "Volcano Plot — Differential Expression",
       x     = "log2 Fold Change",
       y     = "-log10 (Adjusted p-value)",
       color = NULL)

ggsave(file.path(plots_dir, "volcano_plot.png"), p_volcano, width = 10, height = 8)

# ── Plot 3: Heatmap top DE genes ──────────────────────────────────────────────
vsd       <- vst(dds, blind = FALSE)
top_ids   <- rownames(head(res_sig, top_n_heatmap))

if (length(top_ids) > 0) {
  heatmap_mat <- assay(vsd)[top_ids, ]

  sym_hm <- mapIds(org.Hs.eg.db,
                   keys      = top_ids,
                   keytype   = "ENSEMBL",
                   column    = "SYMBOL",
                   multiVals = "first")
  rownames(heatmap_mat) <- make.unique(sym_hm)

  png(file.path(plots_dir, "heatmap_top_DE.png"),
      width = 900, height = 1100, res = 120)
  pheatmap(heatmap_mat,
           scale          = "row",
           annotation_col = metadata,
           main           = paste("Top", top_n_heatmap, "DE genes"),
           color          = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           show_rownames  = TRUE,
           cluster_cols   = TRUE)
  dev.off()
}

cat("\nDE analysis complete. Outputs saved to:", output_dir, "\n")
