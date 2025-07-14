#!/bin/bash


bam_files=(*.duplicates_marked.bam)

process_bam() {
    local bam=$1
    local sample_name=$(basename "$bam" .duplicates_marked.bam)    # Extract the sample name by removing the suffix.
    local count=$(samtools view -c -f1024 "$bam")    # Use samtools to count the number of duplicate reads.
    echo -e "${sample_name}\t${count}"    # Print the sample name and the count, separated by a tab.
}

export -f process_bam

printf "%s\n" "${bam_files[@]}" | parallel -j 20 process_bam
