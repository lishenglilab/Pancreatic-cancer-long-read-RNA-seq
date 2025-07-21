#!/bin/bash

# Set Java memory and Picard path
JAVA_MEM="-Xmx64g"
PICARD_JAR="/public/miniconda3/envs/java_env/share/picard-2.20.4-0/picard.jar"

# Reference file paths
REF_FLAT="/public/reference/star/human/v41/refFlat.txt"
GENOME_FA="/public/reference/star/human/v41/gencode.v41.GRCh38.primary_assembly.genome.fa"

# Find all matching BAM files
BAM_FILES=$(find /public/nanopore_cDNA/bam/ -name "*dedup.sorted.filtered.sorted.bam")

# Check if any files were found
if [ -z "$BAM_FILES" ]; then
    echo "No matching BAM files found"
    exit 1
fi

# Process each BAM file with CollectRnaSeqMetrics
for BAM_FILE in $BAM_FILES; do
    # Extract sample name from BAM filename
    SAMPLE_NAME=$(basename "$BAM_FILE" | cut -d'.' -f1)
    OUTPUT_FILE="${SAMPLE_NAME}_metrics.txt"
    
    echo "Processing: $BAM_FILE"
    echo "Output file: $OUTPUT_FILE"
    
    # Run Picard tool
    java $JAVA_MEM -jar "$PICARD_JAR" \
        CollectRnaSeqMetrics \
        I="$BAM_FILE" \
        O="$OUTPUT_FILE" \
        REF_FLAT="$REF_FLAT" \
        R="$GENOME_FA" \
        STRAND_SPECIFICITY=NONE \
        VALIDATION_STRINGENCY=SILENT \
        METRIC_ACCUMULATION_LEVEL=ALL_READS
    
    if [ $? -eq 0 ]; then
        echo "Successfully processed: $BAM_FILE"
    else
        echo "Failed to process: $BAM_FILE"
    fi
done

echo "All files processed"
