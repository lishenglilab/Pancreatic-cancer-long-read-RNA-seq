#!/bin/bash
mkdir -p /public/nanopore_cDNA/bam/internal_priming

#Parallel Processing
ls /public/nanopore_cDNA/bam/*.bed | parallel -j 20 '
    input_file={}
    base_name=$(basename {} .dedup.sorted.filtered.sorted.bed)
    output_file="/public/nanopore_cDNA/bam/internal_priming/${base_name}_tss.bed"

    awk "BEGIN {OFS=\"\t\"} {
        if (\$6 == \"+\") {
            print \$1, \$2, \$2+1, \$4, \$5, \$6
        } else if (\$6 == \"-\") {
            print \$1, \$3-1, \$3, \$4, \$5, \$6
        }
    }" "$input_file" > "$output_file"

    echo "Processed: $base_name"
'