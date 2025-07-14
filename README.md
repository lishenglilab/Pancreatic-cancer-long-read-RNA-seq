# Analysis of the Nanopore sequencing data pipeline for long RNA in pancreatic cancer

This repository contains a series of scripts to process and analyze ancient DNA (aDNA) sequenced with Oxford Nanopore Technologies. The pipeline starts from raw FASTQ files and performs quality control, contaminant filtering, alignment, post-alignment processing, and downstream analysis including isoform identification and quantification.

The analysis is divided into several major stages, each contained within a numbered directory.

---

## 01: Pre-processing and Quality Control

This stage focuses on assessing the quality of the raw sequencing reads.

- **`0.3.extract_stats.py`**: Parses FASTQ files to calculate read length and average quality score for each read, saving the output to a summary CSV file.
- **`0.1.read_length_plot.py` & `0.2.quality_score_plot.py`**: Python scripts that generate density plots for read length and quality scores from the extracted stats. They help visualize the initial data quality across multiple samples.
- **`0.4.length&quality_density_plot.R`**: An R script that provides an alternative method to generate similar density plots for read length and quality using `ggplot2`.
- **`0.5.NanoReadCategorizer.py`**: Categorizes reads into four bins based on length (<200bp) and quality (<Q7) thresholds, providing a summary report of read characteristics.

## 02: Contaminant Filtering

This stage identifies and removes contaminating reads from the dataset, particularly human DNA.

- **`1.1.kraken.sh`**: Runs Kraken2, a taxonomic classification tool, to identify the species of origin for each read.
- **`1.2.KrakenHumanReadFilter.sh`**: Processes the Kraken2 output to filter out reads classified as human (or other contaminants). It uses `extract_kraken_reads.py` to create new FASTQ files containing only non-human reads.
- **`1.3.kraken_plot.R`**: Generates a stacked bar plot from the Kraken2 results to visualize the taxonomic composition of each sample, showing the proportions of human, bacterial, viral, and other DNA.

## 03: Alignment

The filtered reads are aligned to a reference genome.

- **`2.0.alignment.sh`**: Uses `minimap2` with splicing-aware parameters (`-ax splice`) to align the nanopore reads. It then uses `samtools` to convert, sort, and index the resulting BAM files. The script is designed for use with a SLURM job scheduler.

## 04: Post-alignment Processing

This extensive stage refines the alignment files and performs further quality control.

- **`3.0.bam_filter.sh`**: A simple script to filter out unmapped reads from the BAM files using `samtools view`.
- **Duplicate Marking and Removal**:
    - **`4.0.mark_nanopore_duplicates.py`**: A custom Python script to identify and mark PCR duplicates. Unlike traditional tools, it clusters reads by position and splice junctions, then uses sequence similarity to flag duplicates, which is more suitable for long reads.
    - **`4.1.remove_pcr_duplicates.py`**: Removes the reads that were flagged as duplicates in the previous step.
    - **`4.2.dup_sort_index.sh`**: Sorts and indexes the BAM files that have had duplicates marked.
    - **`4.3.dup_stats.sh`**: Counts the number of marked duplicate reads in each sample.
    - **`4.4.PCR_dep_retio.R`**: Visualizes the duplication rates across all samples with a bar plot.
- **Error Rate Calculation**:
    - **`5.0.compute_nanopore_error_rate.py`**: Calculates the sequencing error rate (mismatches, insertions, deletions) from the alignment information (CIGAR strings and NM tags) in the BAM files.
    - **`5.1.error_rate_calc.sh`**: A wrapper script to run the error rate calculation in parallel for all samples.
    - **`5.2.error_rate_plot.R`**: Generates a bar plot to compare the final error rates across samples.
- **Coverage and Read Counting**:
    - **`6.0.count_reads.sh`**: Counts the total number of reads in the initial FASTQ files.
    - **`6.1.generate_deepth.sh`**: Counts the number of primary mapped reads in the final filtered BAM files.
    - **`6.2.exon_coverage.sh`**: Uses `bedtools` to calculate the per-base coverage over all exonic regions defined in a reference GTF file. This helps assess how well the coding regions of the transcriptome were captured.

## 05: Downstream Analysis (Transcript Isoform Analysis with FLAIR)

This stage focuses on identifying and quantifying transcript isoforms using the FLAIR software.

- **`7.0.bam2bed12.sh`**: Converts the final BAM files into BED12 format, which is required for FLAIR.
- **`7.1.flair_correct.sh`**: Runs `flair correct` to use a reference annotation to correct splice sites in the long reads.
- **`7.2.flair_collapse.sh`**: Runs `flair collapse` to collapse the corrected reads from all samples into a high-confidence set of transcript isoforms.
- **`7.3.flair_quantify_tmp.sh`**: Runs `flair quantify` to calculate the expression levels (in TPM) of the identified isoforms for each sample.

## Other

Miscellaneous analysis scripts.

- **`8.0.pearson.R`**: Reads the TPM expression matrix from FLAIR and generates a Pearson correlation heatmap using `pheatmap`. This is used to visualize the similarity between samples. 
