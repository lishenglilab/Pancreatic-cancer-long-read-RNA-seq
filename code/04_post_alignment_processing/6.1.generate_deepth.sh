#!/bin/bash

# This script counts the number of primary mapped reads in a set of BAM files
# using 'samtools flagstat' and GNU Parallel for concurrent processing.

bam_files=(*.dedup.sorted.filtered.sorted.bam)

process_bam() {
    local bam=$1
    local sample_name=$(basename "$bam" .dedup.sorted.filtered.sorted.bam)

    stats=$(samtools flagstat "$bam")       # Run samtools flagstat to get alignment statistics.

    primary_mapped=$(echo "$stats" | grep "primary mapped")    # Grep the output for the line containing "primary mapped".

    mapped_count=$(echo "$primary_mapped" | awk '{print $1}')    # Use awk to extract the read count (the first field) from that line.

    echo -e "$sample_name\t$mapped_count"
}

export -f process_bam

printf "%s\n" "${bam_files[@]}" | parallel -j 20 process_bam {} > primary_mapped_counts.txt

