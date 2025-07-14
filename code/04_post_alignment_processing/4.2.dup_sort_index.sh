#!/bin/bash


BAM_DIR="/public/nanopore_cDNA/bam"
BAM_FILES=($BAM_DIR/*.duplicates_marked.bam)
INPUT_BAM=${BAM_FILES[$SLURM_ARRAY_TASK_ID-1]}
SAMPLE_NAME=$(basename $INPUT_BAM .duplicates_marked.bam)
SORTED_BAM="$BAM_DIR/${SAMPLE_NAME}.duplicates_marked.sorted.bam"

# Sort the input BAM file. -@ specifies the number of additional threads.
echo "Sorting $INPUT_BAM..."
samtools sort -@ 1 -o $SORTED_BAM $INPUT_BAM

# Create an index (.bai) for the newly sorted BAM file.
echo "Indexing $SORTED_BAM..."
samtools index $SORTED_BAM

echo "Completed processing for sample $SAMPLE_NAME"
