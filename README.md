# Pancreatic-cancer-long-read-RNA-seq

# Nanopore Sequencing Data Analysis Pipeline

This repository contains a series of scripts to process and analyze nanopore sequencing data. The pipeline follows a standard workflow from raw data quality control to downstream analysis.

## Analysis Workflow

The analysis is divided into the following main steps:

1.  [Preprocessing and QC](#1-preprocessing-and-qc)
2.  [Contaminant Filtering](#2-contaminant-filtering)
3.  [Alignment](#3-alignment)
4.  [Post-alignment Processing](#4-post-alignment-processing)
5.  [Downstream Analysis](#5-downstream-analysis)

---

### 1. Preprocessing and QC

**Directory:** `01_preprocessing_and_qc/`

This step focuses on evaluating the quality of the raw nanopore reads.

-   `0.1.read_length_plot.py`: Generates a plot of the read length distribution.
-   `0.2.quality_score_plot.py`: Creates a plot showing the distribution of read quality scores.
-   `0.3.extract_stats.py`: Extracts summary statistics from the `sequencing_summary.txt` file.
-   `0.4.length&quality_density_plot.R`: Generates density plots for read length and quality.
-   `0.5.NanoReadCategorizer.py`: Categorizes reads based on length and quality.

---

### 2. Contaminant Filtering

**Directory:** `02_contaminant_filtering/`

This step identifies and removes contaminating reads from the dataset.

-   `1.1.kraken.sh`: Runs Kraken2 to classify reads and identify contaminants.
-   `1.2.KrakenHumanReadFilter.sh`: Filters out non-human reads based on the Kraken2 classification.
-   `1.3.kraken_plot.R`: Generates plots to visualize the Kraken2 results.

---

### 3. Alignment

**Directory:** `03_alignment/`

This step aligns the filtered reads to a reference genome.

-   `2.0.alignment.sh`: Aligns reads using a mapper like minimap2.
-   `3.0.bam_filter.sh`: Filters the resulting BAM file based on alignment quality or other criteria.

---

### 4. Post-alignment Processing

**Directory:** `04_post_alignment_processing/`

This step processes the alignment files to prepare them for downstream analysis.

-   `4.0.mark_nanopore_duplicates.py` & `4.1.remove_pcr_duplicates.py`: Scripts to identify and remove PCR duplicates, which can be an issue in library preparation.
-   `4.2.dup_sort_index.sh`: Sorts and indexes the BAM file after duplicate removal.
-   `4.3.dup_stats.sh` & `4.4.PCR_dep_retio.R`: Calculate and visualize statistics about the level of PCR duplication.
-   `5.0.compute_nanopore_error_rate.py`, `5.1.error_rate_calc.sh`, & `5.2.error_rate_plot.R`: A set of scripts to calculate and visualize the sequencing error rate from the alignments.

---

### 5. Downstream Analysis

**Directory:** `05_downstream_analysis/`

This final step involves various analyses based on the processed alignments.

-   `6.0.count_reads.sh`: Counts the number of reads in the final BAM file.
-   `6.1.generate_deepth.sh`: Calculates the sequencing depth across the genome.
-   `6.2.exon_coverage.sh`: Calculates coverage over exonic regions.
-   **FLAIR pipeline (`7.*` scripts):**
    -   `7.0.bam2bed12.sh`: Converts BAM to BED12 format.
    -   `7.1.flair_correct.sh`: Uses FLAIR to correct splice sites based on annotations.
    -   `7.2.flair_collapse.sh`: Collapses corrected reads into a set of high-confidence isoforms.
    -   `7.3.flair_quantify_tmp.sh`: Quantifies the expression of the identified isoforms.
-   `8.0.pearson.R`: Performs Pearson correlation analysis, likely comparing sample expressions or other quantitative metrics.

---
## Other

**Directory:** `other/`

This directory may contain miscellaneous scripts or resources not part of the main workflow. 
