#!/bin/bash


bam_files=(*.duplicates_marked.bam)

process_bam() {
    local bam=$1
    local sample_name=$(basename "$bam" .duplicates_marked.bam)
    local count=$(samtools view -c -f1024 "$bam")
    echo -e "${sample_name}\t${count}"
}

export -f process_bam

printf "%s\n" "${bam_files[@]}" | parallel -j 20 process_bam
