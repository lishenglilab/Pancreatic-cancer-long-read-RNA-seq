#!/bin/bash

# This script filters a list of BAM files to remove unmapped reads.
# It reads a list of input BAM file paths from 'bam_list.txt'

cat bam_list.txt | parallel -j 20 "samtools view -@ 3 -b -F 4 {} > {.}.filtered.bam"
