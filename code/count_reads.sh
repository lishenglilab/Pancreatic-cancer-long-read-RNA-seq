#!/bin/bash

INPUT_DIR="/public/nanopore_cDNA/fastq"

FILES=($(find "$INPUT_DIR" -name "filtered_*_nanofilt.fastq"))

parallel -j 20 '
    sample=$(basename {} | sed "s/filtered_\(.*\)_nanofilt.fastq/\1/")
    total_lines=$(wc -l < {} || echo 0)
    reads=$((total_lines / 4))
    echo "$sample,$reads"
' ::: "${FILES[@]}" > reads_count.csv

