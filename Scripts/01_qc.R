# ─────────────────────────────────────────────────────────────────────────────
# 01_qc.R — Quality Control
# RNA-Seq Pipeline | Karim Ghazouani
# ─────────────────────────────────────────────────────────────────────────────
# DESCRIPTION:
#   Runs quality control on raw FASTQ files using fastqcr.
#   Generates per-sample QC reports and a summary MultiQC-style table.
#
# INPUT:
#   data/raw/         — directory containing raw .fastq.gz files
#
# OUTPUT:
#   results/qc/       — FastQC HTML reports + summary table
# ─────────────────────────────────────────────────────────────────────────────

library(fastqcr)

# ── Paths ─────────────────────────────────────────────────────────────────────
raw_dir    <- "data/raw"
output_dir <- "results/qc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Run FastQC on all FASTQ files ─────────────────────────────────────────────
fastq_files <- list.files(raw_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

if (length(fastq_files) == 0) stop("No FASTQ files found in ", raw_dir)

cat("Running FastQC on", length(fastq_files), "files...\n")

fastqc(fq.dir   = raw_dir,
       qc.dir   = output_dir,
       threads  = 4)

# ── Aggregate QC results ──────────────────────────────────────────────────────
qc_summary <- qc_aggregate(output_dir)

cat("\n── QC Summary ──────────────────────────────────────────────────────\n")
print(qc_summary)

# Flag samples with failures
failed <- qc_summary[qc_summary$STATUS == "FAIL", ]
if (nrow(failed) > 0) {
  cat("\n⚠  Samples with QC failures:\n")
  print(failed)
} else {
  cat("\n✓  All samples passed QC.\n")
}

# Save summary to file
write.csv(qc_summary,
          file      = file.path(output_dir, "qc_summary.csv"),
          row.names = FALSE)

cat("\nQC reports saved to:", output_dir, "\n")
