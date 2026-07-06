# ─────────────────────────────────────────────────────────────────────────────
# 03_alignment.R — Read Alignment
# RNA-Seq Pipeline | Karim Ghazouani
# ─────────────────────────────────────────────────────────────────────────────
# DESCRIPTION:
#   Aligns trimmed reads to the reference genome using Rsubread (subread aligner).
#   Builds the genome index if not already present, then aligns all samples.
#
# INPUT:
#   data/trimmed/      — trimmed .fastq.gz files
#   reference/genome.fa — reference genome FASTA
#
# OUTPUT:
#   results/alignment/ — sorted BAM files (.bam) + alignment summary stats
# ─────────────────────────────────────────────────────────────────────────────

library(Rsubread)

# ── Parameters ────────────────────────────────────────────────────────────────
trimmed_dir   <- "data/trimmed"
genome_fasta  <- "reference/genome.fa"
index_dir     <- "reference/subread_index/genome"
output_dir    <- "results/alignment"
threads       <- 4

dir.create(output_dir,         recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(index_dir), recursive = TRUE, showWarnings = FALSE)

# ── Build genome index (only if not already built) ────────────────────────────
index_files <- list.files(dirname(index_dir),
                          pattern = "^genome\\..*\\.b\\.array$")

if (length(index_files) == 0) {
  cat("Building Subread genome index...\n")
  buildindex(basename = index_dir,
             reference = genome_fasta)
  cat("Index built.\n")
} else {
  cat("Genome index already exists. Skipping build.\n")
}

# ── Align each sample ─────────────────────────────────────────────────────────
r1_files     <- sort(list.files(trimmed_dir, pattern = "_R1_trimmed\\.fastq\\.gz$", full.names = TRUE))
r2_files     <- sort(list.files(trimmed_dir, pattern = "_R2_trimmed\\.fastq\\.gz$", full.names = TRUE))
sample_names <- sub("_R1_trimmed.*", "", basename(r1_files))

if (length(r1_files) == 0) stop("No trimmed FASTQ files found in ", trimmed_dir)

paired_end <- length(r2_files) == length(r1_files)
cat("Aligning", length(r1_files), "samples (paired-end:", paired_end, ")\n\n")

stats_list <- list()

for (i in seq_along(r1_files)) {
  name    <- sample_names[i]
  out_bam <- file.path(output_dir, paste0(name, ".bam"))

  cat("Aligning:", name, "\n")

  if (paired_end) {
    stats <- align(index          = index_dir,
                   readfile1      = r1_files[i],
                   readfile2      = r2_files[i],
                   output_file    = out_bam,
                   nthreads       = threads,
                   sortReadsByCoordinates = TRUE)
  } else {
    stats <- align(index          = index_dir,
                   readfile1      = r1_files[i],
                   output_file    = out_bam,
                   nthreads       = threads,
                   sortReadsByCoordinates = TRUE)
  }

  stats_list[[name]] <- stats
  cat("  Done →", out_bam, "\n\n")
}

# ── Save alignment statistics ─────────────────────────────────────────────────
alignment_stats <- do.call(rbind, stats_list)
write.csv(alignment_stats,
          file      = file.path(output_dir, "alignment_stats.csv"),
          row.names = TRUE)

cat("Alignment complete.\n")
cat("BAM files and stats saved to:", output_dir, "\n")
