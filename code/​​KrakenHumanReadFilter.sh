#!/bin/bash

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