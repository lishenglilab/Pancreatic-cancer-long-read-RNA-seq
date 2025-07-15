#!/bin/bash

# This script calculates the per-base read coverage across all exonic regions
# for a set of BAM files. It first generates a BED file of merged exons from a
# GTF annotation, then uses bedtools and GNU Parallel to compute coverage for each sample.

# --- Configuration ---
GTF=/reference/star/human/v41/gencode.v41.primary_assembly.annotation.gtf"
BAM_DIR="/public/nanopore_cDNA/bam"
OUT_DIR="/public/nanopore_cDNA/bam/coverage_results"
mkdir -p $OUT_DIR

# --- Function Definition ---
# This function processes a single BAM file to calculate exon coverage.
process_bam() {
    local bam=$1
    local sample=$(basename "$bam" .dedup.sorted.bam)
    
    # Check for a BAM index (.bai file). If it doesn't exist, create it.
    # This is required by bedtools. It checks for both .bam.bai and .bai extensions.
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        echo "[$(date)] Indexing $bam..."
        if ! samtools index "$bam"; then
            echo "ERROR: Failed to index $bam" >&2
            return 1
        fi
    fi

    # Use bedtools coverage to calculate per-base depth.
    echo "[$(date)] Processing $sample..."
    bedtools coverage -a merged_exons.bed -b "$bam" -d > "$OUT_DIR/${sample}.coverage_per_base.txt" || {
        echo "ERROR: bedtools failed for $bam" >&2
        return 1
    }
}
export -f process_bam
export OUT_DIR

# --- Main Workflow ---
echo "=== Starting Exon Coverage Analysis ==="
date

# Step 1: Create a merged BED file of all exons. This is done only once.
if [[ ! -f "merged_exons.bed" ]]; then
    echo "[$(date)] Extracting exons from GTF..."
    grep -w "exon" "$GTF" | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5, ".", ".", $7}' > exons.bed || exit 1
    bedtools sort -i exons.bed | bedtools merge -s -i - > merged_exons.bed || exit 1
fi

# Step 2: Find and process all BAM files in parallel.
echo "[$(date)] Processing BAM files (20 parallel jobs)..."
find "$BAM_DIR" -name "*.dedup.sorted.bam" | parallel -j 10 --halt soon,fail=1 process_bam || {
    echo "ERROR: Parallel execution failed" >&2
    exit 1
}

