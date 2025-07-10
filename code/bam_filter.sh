#!/bin/bash

cat bam_list.txt | parallel -j 20 "samtools view -@ 3 -b -F 4 {} > {.}.filtered.bam"
