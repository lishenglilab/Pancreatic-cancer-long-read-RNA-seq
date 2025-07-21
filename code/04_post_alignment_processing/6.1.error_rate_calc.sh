#!/bin/bash

# Find all dedup.sorted.bam files (excluding index files)
BAM_FILES=($(ls *.dedup.sorted.bam | grep -v '\.bai$'))

# Run the error rate calculation in parallel
for bam in "${BAM_FILES[@]}"; do
    srun -n1 -c1 --exclusive \
    python /public/nanopore_cDNA/scripts/py/compute_nanopore_error_rate3.py "$bam" \
    > "${bam%.dedup.sorted.bam}_error_rate.txt" &
done

wait
