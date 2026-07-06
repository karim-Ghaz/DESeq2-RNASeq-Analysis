# ─────────────────────────────────────────────────────────────────────────────
# 04_quantification.R — Read Quantification
# RNA-Seq Pipeline | Karim Ghazouani
# ─────────────────────────────────────────────────────────────────────────────
# DESCRIPTION:
#   Counts reads per gene using featureCounts (Rsubread).
#   Produces a raw count matrix ready for normalization and DE analysis.
#
# INPUT:
#   results/alignment/  — sorted BAM files
#   reference/annotation.gtf — gene annotation
#
# OUTPUT:
#   results/counts/counts_matrix.csv — raw gene-level count matrix
#   results/counts/counts_stats.csv  — per-sample counting statistics
# ─────────────────────────────────────────────────────────────────────────────

library(Rsubread)

# ── Parameters ────────────────────────────────────────────────────────────────
alignment_dir  <- "results/alignment"
annotation_gtf <- "reference/annotation.gtf"
output_dir     <- "results/counts"
threads        <- 4
paired_end     <- TRUE   # Set to FALSE for single-end data
strand_specific <- 0     # 0 = unstranded, 1 = stranded, 2 = reversely stranded

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Collect BAM files ─────────────────────────────────────────────────────────
bam_files <- sort(list.files(alignment_dir, pattern = "\\.bam$", full.names = TRUE))

if (length(bam_files) == 0) stop("No BAM files found in ", alignment_dir)

cat("Quantifying", length(bam_files), "samples...\n")

# ── Run featureCounts ─────────────────────────────────────────────────────────
fc <- featureCounts(
  files            = bam_files,
  annot.ext        = annotation_gtf,
  isGTFAnnotationFile = TRUE,
  GTF.featureType  = "exon",
  GTF.attrType     = "gene_id",
  isPairedEnd      = paired_end,
  strandSpecific   = strand_specific,
  nthreads         = threads,
  countMultiMappingReads = FALSE,
  requireBothEndsMapped  = paired_end
)

# ── Extract and clean count matrix ───────────────────────────────────────────
counts_matrix <- fc$counts

# Clean column names: remove path and .bam extension
colnames(counts_matrix) <- sub("\\.bam$", "",
                                basename(colnames(counts_matrix)))

cat("\nCount matrix dimensions:", nrow(counts_matrix), "genes x",
    ncol(counts_matrix), "samples\n")

# ── Save outputs ──────────────────────────────────────────────────────────────
write.csv(counts_matrix,
          file      = file.path(output_dir, "counts_matrix.csv"),
          row.names = TRUE)

write.csv(fc$stat,
          file      = file.path(output_dir, "counts_stats.csv"),
          row.names = FALSE)

cat("\nCount matrix saved to:", file.path(output_dir, "counts_matrix.csv"), "\n")
cat("Counting stats saved to:", file.path(output_dir, "counts_stats.csv"), "\n")
