#!/bin/bash


manifest="/public/nanopore_cDNA/flair/flair/reads_manifest.tsv"
isoforms="/public/nanopore_cDNA/flair/flair/all_samples_collapsed3.isoforms.fa"
out_dir="/public/nanopore_cDNA/flair/quantify/"
out_name=${out_dir}"flair_quantify"

# Run flair quantify:
flair quantify -r ${manifest} -i ${isoforms} --output ${out_name} --threads 72 --temp_dir ${out_dir} --sample_id_only --tpm

