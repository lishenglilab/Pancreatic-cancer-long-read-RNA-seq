#!/bin/bash

# This script filters FASTQ files based on Kraken2 classification results.
# It uses 'extract_kraken_reads.py' to select reads classified as Human (taxid 9606)
# or unclassified (taxid 0), and saves them to new FASTQ files.
# The script leverages GNU Parallel to process multiple samples simultaneously.

# --- Configuration ---
FASTQ_DIR="/data/transcript/fq"
KRAKEN_DIR="kraken_standard_out"
OUTPUT_DIR="filtered_fastq_human"
NUM_JOBS=20

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# --- Main Logic ---

# Function to process a single sample
process_sample() {
    output_file="$1"
    sample=$(basename "$output_file" _output.txt)
    fastq_file="$FASTQ_DIR/$sample/$sample.clean.fastq.gz"
    filtered_file="$OUTPUT_DIR/$sample/filtered_$sample.fastq"

    if [ -f "$fastq_file" ]; then
        echo "Processing sample: $sample"
        extract_kraken_reads.py -k "$output_file" \
                                -s "$fastq_file" \
                                -t 9606 0 \
                                --fastq-output \
                                -o "$filtered_file"
    else
        echo "Error: $fastq_file does not exist, skipping sample $sample"
    fi
}
export -f process_sample

# Use GNU Parallel to process all Kraken2 output files in parallel
find "$KRAKEN_DIR" -name "*_output.txt" | parallel -j 20 process_sample {}

echo "All samples have been processed! Filtered files have been saved to $OUTPUT_DIR"
