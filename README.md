# Analysis of the Nanopore Sequencing Data Pipeline for Long RNA in Pancreatic Cancer

This repository contains a comprehensive pipeline for processing and analyzing long-read RNA sequencing data from Oxford Nanopore Technologies, specifically focused on pancreatic cancer samples. The pipeline covers the entire workflow from raw FASTQ files through quality control, contaminant filtering, alignment, post-alignment processing, and downstream analysis including transcript isoform identification and quantification.

## Pipeline Overview

The analysis is organized into several sequential stages, each contained within a numbered directory:

1. Pre-processing and Quality Control
2. Filtering
3. Alignment to Reference Genome
4. Post-alignment Processing
5. Downstream Analysis (Transcript Isoform Analysis)

## Software Requirements

This pipeline utilizes several bioinformatics tools. The specific versions used are documented in `code/version.txt` and include:

- Kraken2 (v2.1.3) - Taxonomic classification
- NanoFilt (v2.8.0) - Read filtering by length and quality
- minimap2 (v2.22-r1101) - Long-read alignment
- samtools (v1.13) - SAM/BAM file manipulation
- FLAIR (v2.0.0) - Full-Length Alternative Isoform analysis of RNA
- picard(2.20.4-SNAPSHOT) - (CollectRnaSeqMetrics) calculates RNA sequencing quality metrics

## Detailed Pipeline Description

### 01: Pre-processing and Quality Control

This stage focuses on assessing the quality of the sequencing reads.

- **`1.0.extract_stats.py`**: Parses FASTQ files to calculate read length and average quality score for each read, saving the output to a summary CSV file. Uses a mathematically correct method for averaging Phred quality scores by converting to error probabilities.
- **`1.1.length&quality_density_plot.R`**: Generates density plots for read length and quality distributions using `ggplot2`.

### 02: Filtering

This stage identifies and removes contaminating reads from the dataset.

- **`2.0.run_kraken.sh`**: Runs Kraken2 taxonomic classification tool to identify the species of origin for each read.
- **`2.1.KrakenHumanReadFilter.sh`**: Process the Kraken2 output results, retaining sequences identified as "Human" or "Unclassified," and generate a new FASTQ file containing only the target reads.
- **`2.2.kraken_plot.R`**: Generates a stacked bar plot from the Kraken2 results to visualize the taxonomic composition of each sample.
- **`2.3.NanoReadCategorizer.py`**: Categorizes reads into four bins based on length (<200bp) and quality (<Q7) thresholds, providing a summary report of read characteristics.
- **`2.4.filter_reads_by_length_quality.sh`**: Filters reads based on minimum length (200bp) and quality score (Q7) thresholds using NanoFilt.

### 03: Alignment

The filtered reads are aligned to a reference genome.

- **`3.0.alignment.sh`**: Uses `minimap2` with splicing-aware parameters (`-ax splice`) to align the nanopore reads to a reference genome. The script then uses `samtools` to convert, sort, and index the resulting BAM files. Designed for use with a SLURM job scheduler.

### 04: Post-alignment Processing

This extensive stage refines the alignment files and performs further quality control.

- **`4.0.bam_filter.sh`**: Filters out unmapped reads from the BAM files using `samtools view`.
- **Duplicate Marking and Removal**:
    - **`5.0.mark_nanopore_duplicates.py`**: Used to identify and mark PCR duplicate sequences. This script clusters reads by position and splice sites, then marks duplicates based on sequence similarity.
    - **`5.1.remove_pcr_duplicates.py`**: Removes the reads that were flagged as duplicates in the previous step.
    - **`5.2.dup_sort_index.sh`**: Sorts and indexes the BAM files after duplicate removal.
    - **`5.3.dup_stats.sh`**: Counts the number of marked duplicate reads in each sample.
    - **`5.4.PCR_dep_retio.R`**: Visualizes the duplication rates across all samples with a bar plot.
- **Error Rate Calculation**:
    - **`6.0.compute_nanopore_error_rate.py`**: Calculates the sequencing error rate (mismatches, insertions, deletions) from the alignment information (CIGAR strings and NM tags) in the BAM files.
    - **`6.1.error_rate_calc.sh`**: A wrapper script to run the error rate calculation in parallel for all samples.
    - **`6.2.error_rate_plot.R`**: Generates a bar plot to compare the final error rates across samples.
- **Coverage and Read Counting**:
    - **`7.0.count_reads.sh`**: Counts the total number of reads in the initial FASTQ files.
    - **`7.1.generate_deepth.sh`**: Counts the number of primary mapped reads in the final filtered BAM files.
    - **`7.2.exon_coverage.sh`**: Automatically process multiple RNA-seq BAM files using Picard's CollectRnaSeqMetrics tool to generate quality metric reports for each sample.
    - **`7.3.plot_gene_body_coverage.R`**: Visualizes the coverage distribution across gene bodies.
- **Internal Priming Analysis**:
    - **`8.0.extract_tss_from_bed.sh`**: Extracts transcription start sites (TSS) from BED files.
    - **`8.1.Genomic_seq_star.sh`**: Retrieves genomic sequences around the TSS.
    - **`8.2.internal_priming_analysis.sh`**: Analyzes sequences for internal priming artifacts.

### 05: Downstream Analysis (Transcript Isoform Analysis with FLAIR)

This stage focuses on identifying and quantifying transcript isoforms using the FLAIR software.

- **`9.0.bam2bed12.sh`**: Converts the final BAM files into BED12 format, which is required for FLAIR.
- **`9.1.flair_correct.sh`**: Runs `flair correct` to use a reference annotation to correct splice sites in the long reads.
- **`9.2.flair_collapse.sh`**: Runs `flair collapse` to collapse the corrected reads from all samples into a high-confidence set of transcript isoforms.
- **`9.3.flair_quantify_tmp.sh`**: Runs `flair quantify` to calculate the expression levels (in TPM) of the identified isoforms for each sample.

### Other

Miscellaneous analysis scripts.

- **`10.0.pearson.R`**: Reads the TPM expression matrix from FLAIR and generates a Pearson correlation heatmap using `pheatmap`. This visualizes the similarity between samples.
