# ─────────────────────────────────────────────────────────────────────────────
# 07_functional_analysis.R — Functional Enrichment Analysis
# RNA-Seq Pipeline | Dr. Hicham Charoute
# ─────────────────────────────────────────────────────────────────────────────
# DESCRIPTION:
#   Performs Gene Ontology (GO) and KEGG pathway enrichment analysis on
#   differentially expressed genes using clusterProfiler.
#   Produces dotplot, barplot, and enrichment map visualizations.
#
# INPUT:
#   results/de/de_results_significant.csv — significant DE genes
#
# OUTPUT:
#   results/functional/go_bp_results.csv      — GO Biological Process results
#   results/functional/go_mf_results.csv      — GO Molecular Function results
#   results/functional/kegg_results.csv       — KEGG pathway results
#   results/functional/plots/                 — enrichment plots
# ─────────────────────────────────────────────────────────────────────────────

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# ── Parameters ────────────────────────────────────────────────────────────────
de_results_file <- "results/de/de_results_significant.csv"
output_dir      <- "results/functional"
plots_dir       <- file.path(output_dir, "plots")
padj_method     <- "BH"
show_categories <- 20

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load DE results ───────────────────────────────────────────────────────────
de_results  <- read.csv(de_results_file, row.names = 1)
ensembl_ids <- rownames(de_results)

if (length(ensembl_ids) == 0) stop("No significant DE genes found.")
cat("Input:", length(ensembl_ids), "significant DE genes\n")

# ── Convert Ensembl IDs to Entrez IDs ─────────────────────────────────────────
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys    = ensembl_ids,
                     column  = "ENTREZID",
                     keytype = "ENSEMBL")

entrez_ids <- na.omit(unique(entrez_ids))
cat("Mapped to Entrez:", length(entrez_ids), "IDs\n\n")

# ── GO Enrichment — Biological Process ────────────────────────────────────────
cat("Running GO enrichment (Biological Process)...\n")
go_bp <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = padj_method,
                  readable      = TRUE)

# ── GO Enrichment — Molecular Function ────────────────────────────────────────
cat("Running GO enrichment (Molecular Function)...\n")
go_mf <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = padj_method,
                  readable      = TRUE)

# ── KEGG Pathway Enrichment ───────────────────────────────────────────────────
cat("Running KEGG pathway enrichment...\n")
kegg <- enrichKEGG(gene          = entrez_ids,
                   organism      = "hsa",     # hsa = Homo sapiens
                   pAdjustMethod = padj_method)

# ── Save results ──────────────────────────────────────────────────────────────
if (!is.null(go_bp) && nrow(go_bp) > 0) {
  write.csv(as.data.frame(go_bp),
            file.path(output_dir, "go_bp_results.csv"), row.names = FALSE)
  cat("GO-BP: found", nrow(go_bp), "enriched terms\n")
} else { cat("GO-BP: no significant terms\n") }

if (!is.null(go_mf) && nrow(go_mf) > 0) {
  write.csv(as.data.frame(go_mf),
            file.path(output_dir, "go_mf_results.csv"), row.names = FALSE)
  cat("GO-MF: found", nrow(go_mf), "enriched terms\n")
} else { cat("GO-MF: no significant terms\n") }

if (!is.null(kegg) && nrow(kegg) > 0) {
  write.csv(as.data.frame(kegg),
            file.path(output_dir, "kegg_results.csv"), row.names = FALSE)
  cat("KEGG: found", nrow(kegg), "enriched pathways\n")
} else { cat("KEGG: no significant pathways\n") }

# ── Plots ─────────────────────────────────────────────────────────────────────

# GO-BP dotplot
if (!is.null(go_bp) && nrow(go_bp) > 0) {
  p1 <- dotplot(go_bp, showCategory = show_categories) +
    ggtitle("GO Enrichment — Biological Process")
  ggsave(file.path(plots_dir, "go_bp_dotplot.png"), p1, width = 10, height = 9)

  p2 <- barplot(go_bp, showCategory = show_categories) +
    ggtitle("GO Enrichment — Biological Process")
  ggsave(file.path(plots_dir, "go_bp_barplot.png"), p2, width = 10, height = 9)
}

# GO-MF dotplot
if (!is.null(go_mf) && nrow(go_mf) > 0) {
  p3 <- dotplot(go_mf, showCategory = show_categories) +
    ggtitle("GO Enrichment — Molecular Function")
  ggsave(file.path(plots_dir, "go_mf_dotplot.png"), p3, width = 10, height = 8)
}

# KEGG dotplot
if (!is.null(kegg) && nrow(kegg) > 0) {
  p4 <- dotplot(kegg, showCategory = show_categories) +
    ggtitle("KEGG Pathway Enrichment")
  ggsave(file.path(plots_dir, "kegg_dotplot.png"), p4, width = 10, height = 8)
}

cat("\nFunctional analysis complete. Outputs saved to:", output_dir, "\n")
