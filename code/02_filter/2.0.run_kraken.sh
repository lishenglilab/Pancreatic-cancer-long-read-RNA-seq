#!/bin/bash

# This script runs the Kraken2 tool for taxonomic classification on multiple FASTQ files.
# It iterates through sample directories, runs Kraken2 on gzipped FASTQ files,
# and saves the output and report for each sample.


DATA_DIR="/data/transcript/fq"
KRAKEN_DB="k2_standard_20241228"
OUTPUT_DIR="kraken_standard_out"
THREADS=20

# Create output directory if it doesn't exist
mkdir -p kraken_standard_out

# Iterate through all sample directories
for dir in /data/yvzeng/Project/pancreasLongRead/data/transcript/fq/*/; do
    sample_name=$(basename "${dir%/}")

    # Run Kraken2 analysis
    kraken2 --db k2_standard_20241228 \
            --threads 20 \
            --output "kraken_standard_out/${sample_name}_output.txt" \
            --report "kraken_standard_out/${sample_name}_report.txt" \
            "${dir}clean.fastq.gz"
done
