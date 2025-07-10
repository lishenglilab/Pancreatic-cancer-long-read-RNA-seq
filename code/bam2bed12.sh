#!/bin/bash


find . -name "*.dedup.sorted.filtered.sorted.bam" | parallel -j 20 "samtools index {}"

find . -name "*.dedup.sorted.filtered.sorted.bam" | while read bam; do
    echo "python /public/software/flair-2.0.0/src/flair/bam2Bed12.py -i \"$bam\" --keep_supplementary > \"${bam%.bam}.bed\""
done > joblist.txt

parallel -j $SLURM_CPUS_PER_TASK < joblist.txt
