#!/bin/bash


bam_files=(*.dedup.sorted.filtered.sorted.bam)

process_bam() {
    local bam=$1
    local sample_name=$(basename "$bam" .dedup.sorted.filtered.sorted.bam)

    stats=$(samtools flagstat "$bam")

    primary_mapped=$(echo "$stats" | grep "primary mapped")

    mapped_count=$(echo "$primary_mapped" | awk '{print $1}')

    echo -e "$sample_name\t$mapped_count"
}

export -f process_bam

printf "%s\n" "${bam_files[@]}" | parallel -j 20 process_bam {} > primary_mapped_counts.txt

