#!/bin/bash
cd /public/nanopore_cDNA/bam/internal_priming


# Process all *_tss.bed files in parallel (20 jobs at a time)
ls *_tss.bed | parallel -j 20 '

sample=$(basename {} _tss.bed)
# Add 20bp flanking regions to each side of the TSS coordinates
bedtools slop -i {} -g genome.chrom.sizes -b 20 | \
# Extract the DNA sequences for these regions from reference genome
bedtools getfasta -fi /public/home/shengli-lab/shengnanluo/reference/star/human/v41/gencode.v41.GRCh38.primary_assembly.genome.fa -bed - > ${sample}_flanks.fa
echo "Processing completed.: {}"
