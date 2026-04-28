# RNA-Seq Pipeline: Zika Virus Infection in hNPCs

This repository contains a full end-to-end RNA-sequencing data analysis pipeline developed for the BioProject **PRJNA313294**. 

## Project Objective
The primary goal of this project is to evaluate the technical concordance between two sequencing platforms: **MiSeq** (Paired-End 75nt) and **NextSeq** (Single-End 75nt). The analysis checks whether the biological signal (Zika infection vs. Mock control) remains robust and distinct regardless of the sequencing device used.

## Dataset Details
* **Organism/Cell Type:** Human Neural Progenitor Cells (hNPCs)
* **Conditions:** Zika-infected (2 samples) vs. Mock control (2 samples)
* **Libraries:**
  * 4x MiSeq libraries (75 nt Paired-End)
  * 4x NextSeq libraries (75 nt Single-End)
* **Reference Genome:** Human `hg19`

## Workflow Structure

### 1. Upstream Processing (Bash Pipeline)
The `pipeline.sh` script automates the retrieval and processing of raw sequencing data:
* **Data Retrieval:** Downloads SRR files using `fastq-dump`.
* **Quality Control:** Runs `FastQC`.
* **Trimming:** Uses `Trimmomatic` to remove adapters and low-quality bases.
* **Mapping:** Aligns reads to the `hg19` reference genome using `HISAT2`.
* **Format Conversion:** Converts, sorts, and indexes SAM to BAM files via `SAMtools`.
* **Quantification:** Generates feature counts using `FeatureCounts`.

### 2. Downstream Analysis (RMarkdown)
The `.Rmd` script utilizes `DESeq2` to analyze the generated count matrices (`counts_SE.txt` and `counts_PE.txt`):
* **Data merging & normalization** (rlog transformation).
* **Correlation Analysis:** Pearson correlation confirms highly consistent expression levels between MiSeq and NextSeq (r = 0.99).
* **Visualization:**
  * Heatmaps (`pheatmap`, `gplots`) showing the top 1000 highly variable genes and highly expressed genes.
  * PCA plot demonstrating that samples cluster strictly by biological condition rather than the sequencing platform.
