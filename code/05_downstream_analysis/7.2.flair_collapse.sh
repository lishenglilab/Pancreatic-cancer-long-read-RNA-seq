#!/bin/bash

# This script runs the 'flair collapse' step of the flair pipeline.

# --- Configuration ---
ref_gtf="/reference/star/human/v41/gencode.v41.primary_assembly.annotation.gtf"
ref_fa="/reference/star/human/v41/gencode.v41.GRCh38.primary_assembly.genome.fa"
bed_file="/public/nanopore_cDNA/filar/ffilar/merged_corrected.bed"
output_dir="/public/nanopore_cDNA/filar/ffilar"
fastq_dir="/public/nanopore_cDNA/fastq"

# Create output directories
mkdir -p "${output_dir}"
mkdir -p "${output_dir}/tmp3"

# Find all the original filtered FASTQ files.
fastq_files=()
for sample in $(ls "${fastq_dir}"); do
    fastq_files+=("${fastq_dir}/${sample}/filtered_${sample}_nanofilt.fastq")
done

fastq_list=$(printf "%s," "${fastq_files[@]}" | sed 's/,$//')



# --- Run flair collapse ---
flair collapse -g "${ref_fa}" -q "${bed_file}" --gtf "${ref_gtf}" \
     -r "${fastq_list}" \
     --output "${output_dir}/all_samples_collapsed3" \
     --threads 50 \
     --temp_dir "${output_dir}/tmp3"

echo "Flair collapse step completed."
