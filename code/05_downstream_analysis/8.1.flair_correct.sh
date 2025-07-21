#!/bin/bash

# This script runs the 'flair correct' step on a set of BED12 files.
# This step uses splice site information from a reference GTF file to
# correct inaccuracies in the splice junctions of the nanopore reads.

# --- Configuration ---
ref_gtf="/reference/star/human/v41/gencode.v41.primary_assembly.annotation.gtf"
ref_fa="/reference/star/human/v41/gencode.v41.GRCh38.primary_assembly.genome.fa"
output_dir="/public/nanopore_cDNA/filar"

# --- Main Logic ---
bed_files=(*dedup.sorted.filtered.sorted.bed)

process_bed() {
    local bed_file=$1
    local sample_name=$(basename "$bed_file" .dedup.sorted.filtered.sorted.bed)
    local out_name="${output_dir}/${sample_name}_corrected"
    
    echo "Processing $bed_file..."
    flair correct -q "$bed_file" -f "$ref_gtf" -g "$ref_fa" --output "$out_name" --threads 3        # Run flair correct:
    echo "Finished processing $bed_file"
}

# Export variables and the function to make them available to GNU Parallel.
export ref_gtf ref_fa output_dir
export -f process_bed

# Use parallel to run the process_bed function on each file.
parallel -j 20 process_bed ::: "${bed_files[@]}"

echo "All flair correct jobs completed"
