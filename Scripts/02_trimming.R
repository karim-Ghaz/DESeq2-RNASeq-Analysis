# ─────────────────────────────────────────────────────────────────────────────
# 02_trimming.R — Read Trimming and Filtering
# RNA-Seq Pipeline | Dr. Hicham Charoute
# ─────────────────────────────────────────────────────────────────────────────
# DESCRIPTION:
#   Trims adapter sequences and low-quality bases from raw reads using
#   Trimmomatic via system calls. Supports both single-end and paired-end data.
#
# INPUT:
#   data/raw/          — raw .fastq.gz files
#   data/adapters.fa   — adapter sequences FASTA (Illumina TruSeq by default)
#
# OUTPUT:
#   data/trimmed/      — trimmed .fastq.gz files
# ─────────────────────────────────────────────────────────────────────────────

# ── Parameters ────────────────────────────────────────────────────────────────
raw_dir      <- "data/raw"
trimmed_dir  <- "data/trimmed"
adapter_file <- "data/adapters.fa"   # Path to adapter FASTA
threads      <- 4

# Trimmomatic parameters
leading      <- 3    # Remove leading low-quality bases below quality 3
trailing     <- 3    # Remove trailing low-quality bases below quality 3
sliding_win  <- "4:15"  # Sliding window: window size 4, quality threshold 15
min_len      <- 36   # Minimum read length after trimming

dir.create(trimmed_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper: run Trimmomatic for one sample ────────────────────────────────────
trim_sample <- function(r1, r2 = NULL, sample_name) {

  if (!is.null(r2)) {
    # Paired-end
    r1_out  <- file.path(trimmed_dir, paste0(sample_name, "_R1_trimmed.fastq.gz"))
    r2_out  <- file.path(trimmed_dir, paste0(sample_name, "_R2_trimmed.fastq.gz"))
    r1_unp  <- file.path(trimmed_dir, paste0(sample_name, "_R1_unpaired.fastq.gz"))
    r2_unp  <- file.path(trimmed_dir, paste0(sample_name, "_R2_unpaired.fastq.gz"))

    cmd <- paste(
      "trimmomatic PE -threads", threads,
      r1, r2,
      r1_out, r1_unp, r2_out, r2_unp,
      paste0("ILLUMINACLIP:", adapter_file, ":2:30:10"),
      paste0("LEADING:", leading),
      paste0("TRAILING:", trailing),
      paste0("SLIDINGWINDOW:", sliding_win),
      paste0("MINLEN:", min_len)
    )
  } else {
    # Single-end
    out <- file.path(trimmed_dir, paste0(sample_name, "_trimmed.fastq.gz"))

    cmd <- paste(
      "trimmomatic SE -threads", threads,
      r1, out,
      paste0("ILLUMINACLIP:", adapter_file, ":2:30:10"),
      paste0("LEADING:", leading),
      paste0("TRAILING:", trailing),
      paste0("SLIDINGWINDOW:", sliding_win),
      paste0("MINLEN:", min_len)
    )
  }

  cat("Trimming:", sample_name, "\n")
  ret <- system(cmd)
  if (ret != 0) warning("Trimmomatic failed for sample: ", sample_name)
}

# ── Detect paired-end or single-end samples ───────────────────────────────────
r1_files <- sort(list.files(raw_dir, pattern = "_R1.*\\.fastq\\.gz$", full.names = TRUE))
r2_files <- sort(list.files(raw_dir, pattern = "_R2.*\\.fastq\\.gz$", full.names = TRUE))

if (length(r1_files) == 0) stop("No FASTQ files found in ", raw_dir)

if (length(r2_files) == length(r1_files)) {
  cat("Detected paired-end data (", length(r1_files), "samples)\n")
  sample_names <- sub("_R1.*", "", basename(r1_files))
  mapply(trim_sample, r1_files, r2_files, sample_names)
} else {
  cat("Detected single-end data (", length(r1_files), "samples)\n")
  sample_names <- sub("_R1.*|\\.fastq\\.gz", "", basename(r1_files))
  mapply(trim_sample, r1_files, MoreArgs = list(r2 = NULL), sample_names)
}

cat("\nTrimming complete. Files saved to:", trimmed_dir, "\n")
