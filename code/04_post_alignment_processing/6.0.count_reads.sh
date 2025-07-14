#!/bin/bash

# This script counts the number of reads in multiple FASTQ files and outputs the
# results to a CSV file. It leverages GNU Parallel to process files concurrently.
# The script assumes a standard FASTQ format where each read occupies 4 lines.

INPUT_DIR="/public/nanopore_cDNA/fastq"

FILES=($(find "$INPUT_DIR" -name "filtered_*_nanofilt.fastq"))

parallel -j 20 '
    sample=$(basename {} | sed "s/filtered_\(.*\)_nanofilt.fastq/\1/")
    total_lines=$(wc -l < {} || echo 0)
    reads=$((total_lines / 4))    # Calculate the number of reads (assuming 4 lines per read).
    echo "$sample,$reads"
' ::: "${FILES[@]}" > reads_count.csv

