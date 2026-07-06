# RNA-Seq Pipeline in R

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.17-green)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Welcome to the **RNA-Seq Pipeline**! This pipeline is designed to take you from raw sequencing data all the way to biological insights using R. It covers all essential steps — from quality control to functional enrichment — ensuring reproducibility and ease of use for your RNA-Seq projects.

---

## Table of Contents

1. [Biological Overview](#biological-overview)
2. [Pipeline Steps](#pipeline-steps)
   - [1. Quality Control](#1-quality-control)
   - [2. Read Trimming and Filtering](#2-read-trimming-and-filtering)
   - [3. Alignment](#3-alignment)
   - [4. Quantification](#4-quantification)
   - [5. Normalization](#5-normalization)
   - [6. Differential Expression Analysis](#6-differential-expression-analysis)
   - [7. Functional Analysis](#7-functional-analysis)
3. [Repository Structure](#repository-structure)
4. [Usage Instructions](#usage-instructions)
5. [Dependencies](#dependencies)
6. [License](#license)
7. [Contact](#contact)

---

## Biological Overview

RNA sequencing (RNA-Seq) is a powerful technique that captures a snapshot of the entire transcriptome — the complete set of RNA transcripts — in a cell at a given moment. By analyzing RNA-Seq data, we can uncover gene expression patterns, identify differentially expressed genes, and understand the molecular mechanisms underlying biological processes and diseases.

This pipeline guides you through each critical step of RNA-Seq data analysis, ensuring accurate and biologically meaningful results.

---

## Pipeline Steps

### 1. Quality Control

**Biological meaning:**
Before any analysis, it is essential to verify the quality of raw sequencing data. Poor-quality reads can introduce noise and bias downstream results.

**What it does:**
Runs FastQC on all raw FASTQ files and aggregates per-sample metrics into a summary report. Flags samples with adapter contamination, low Q-scores, or GC-content anomalies.

**Script:** `scripts/01_qc.R`

---

### 2. Read Trimming and Filtering

**Biological meaning:**
Adapter sequences and low-quality bases do not originate from the biological sample — they are technical artifacts that reduce mapping accuracy and inflate noise.

**What it does:**
Calls Trimmomatic to remove adapters and trim low-quality bases. Supports both paired-end and single-end libraries. Discards reads shorter than a minimum length threshold.

**Script:** `scripts/02_trimming.R`

---

### 3. Alignment

**Biological meaning:**
Mapping reads to a reference genome assigns each read to its genomic origin, which is the foundation for measuring gene expression.

**What it does:**
Builds a Subread genome index (if not already present) and aligns trimmed reads using `Rsubread::align()`. Outputs coordinate-sorted BAM files and per-sample alignment statistics.

**Script:** `scripts/03_alignment.R`

---

### 4. Quantification

**Biological meaning:**
Counting reads overlapping each gene converts alignment data into a numerical matrix of expression values — one value per gene per sample.

**What it does:**
Runs `featureCounts` on all BAM files using the provided GTF annotation. Produces a raw count matrix and per-sample counting statistics.

**Script:** `scripts/04_quantification.R`

---

### 5. Normalization

**Biological meaning:**
Raw counts are affected by technical factors such as sequencing depth and library size. Normalization corrects for these to allow fair comparison between samples.

**What it does:**
Uses DESeq2 median-of-ratios normalization and variance-stabilizing transformation (VST). Produces boxplots and a PCA plot to verify that normalization removes technical variation.

**Script:** `scripts/05_normalization.R`

---

### 6. Differential Expression Analysis

**Biological meaning:**
Identifying genes that are significantly up- or downregulated between conditions reveals the molecular basis of the biological response under study.

**What it does:**
Runs the full DESeq2 pipeline (Wald test). Annotates results with HGNC gene symbols. Produces an MA plot, volcano plot, and heatmap of the top 50 differentially expressed genes.

**Script:** `scripts/06_differential_expression.R`

---

### 7. Functional Analysis

**Biological meaning:**
Individual gene lists gain biological interpretation when placed in the context of pathways and cellular functions. Enrichment analysis reveals which biological processes are most affected.

**What it does:**
Runs GO enrichment (Biological Process and Molecular Function) and KEGG pathway analysis using `clusterProfiler`. Produces dotplots, barplots, and saves enrichment tables.

**Script:** `scripts/07_functional_analysis.R`

---

## Repository Structure

```
RNAseq-Pipeline/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── scripts/
│   ├── 01_qc.R
│   ├── 02_trimming.R
│   ├── 03_alignment.R
│   ├── 04_quantification.R
│   ├── 05_normalization.R
│   ├── 06_differential_expression.R
│   └── 07_functional_analysis.R
│
├── data/
│   ├── raw/            ← place your raw .fastq.gz files here
│   ├── trimmed/        ← generated by 02_trimming.R
│   ├── samples.txt     ← sample metadata (sample | condition)
│   └── adapters.fa     ← adapter sequences for trimming
│
├── reference/
│   ├── genome.fa       ← reference genome FASTA
│   └── annotation.gtf  ← gene annotation GTF
│
└── results/
    ├── qc/
    ├── alignment/
    ├── counts/
    ├── normalization/
    ├── de/
    └── functional/
```

> **Note:** `data/raw/`, `reference/`, and `results/` are not versioned (see `.gitignore`). Add your own files locally.

---

## Usage Instructions

**1. Clone the repository:**

```bash
git clone https://github.com/<your-username>/RNAseq-Pipeline.git
cd RNAseq-Pipeline
```

**2. Prepare your data:**

- Place raw `.fastq.gz` files in `data/raw/`
- Place `genome.fa` and `annotation.gtf` in `reference/`
- Fill in `data/samples.txt` (tab-separated: `sample` and `condition` columns)

**3. Install R dependencies:**

```r
install.packages(c("ggplot2", "ggrepel", "pheatmap",
                   "RColorBrewer", "reshape2", "dplyr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "Rsubread", "fastqcr",
                       "clusterProfiler", "enrichplot",
                       "org.Hs.eg.db"))
```

**4. Run the pipeline step by step:**

```bash
Rscript scripts/01_qc.R
Rscript scripts/02_trimming.R
Rscript scripts/03_alignment.R
Rscript scripts/04_quantification.R
Rscript scripts/05_normalization.R
Rscript scripts/06_differential_expression.R
Rscript scripts/07_functional_analysis.R
```

**5. Review results:**

| Step | Output location |
|---|---|
| QC reports | `results/qc/` |
| BAM files | `results/alignment/` |
| Count matrix | `results/counts/` |
| Normalized counts + PCA | `results/normalization/` |
| DE results + plots | `results/de/` |
| Enrichment tables + plots | `results/functional/` |

---

## Dependencies

**R ≥ 4.2**

| Package | Source | Purpose |
|---|---|---|
| fastqcr | CRAN | Quality control |
| Rsubread | Bioconductor | Alignment + quantification |
| DESeq2 | Bioconductor | Normalization + DE analysis |
| clusterProfiler | Bioconductor | GO and KEGG enrichment |
| enrichplot | Bioconductor | Enrichment visualization |
| org.Hs.eg.db | Bioconductor | Human gene annotation |
| ggplot2 | CRAN | Visualization |
| ggrepel | CRAN | Volcano plot labels |
| pheatmap | CRAN | Heatmap |
| RColorBrewer | CRAN | Color palettes |
| reshape2 | CRAN | Data reshaping |
| dplyr | CRAN | Data manipulation |

**External tools (must be in PATH):**
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Contact

**Karim Ghazouani**

For questions or contributions, feel free to open an issue or submit a pull request.

---

*Keep experimenting, and feel free to adapt the scripts to your own experimental design!* 🧬
