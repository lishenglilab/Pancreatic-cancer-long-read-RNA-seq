#!/bin/bash


ref_gtf="/public/home/shengli-lab/shengnanluo/reference/star/human/v41/gencode.v41.primary_assembly.annotation.gtf"
ref_fa="/public/home/shengli-lab/shengnanluo/reference/star/human/v41/gencode.v41.GRCh38.primary_assembly.genome.fa"
output_dir="/public/nanopore_cDNA/filar"

bed_files=(*dedup.sorted.filtered.sorted.bed)

process_bed() {
    local bed_file=$1
    local sample_name=$(basename "$bed_file" .dedup.sorted.filtered.sorted.bed)
    local out_name="${output_dir}/${sample_name}_corrected"
    
    echo "Processing $bed_file..."
    flair correct -q "$bed_file" -f "$ref_gtf" -g "$ref_fa" --output "$out_name" --threads 3
    echo "Finished processing $bed_file"
}

export ref_gtf ref_fa output_dir
export -f process_bed

parallel -j 20 process_bed ::: "${bed_files[@]}"

echo "All flair correct jobs completed"
