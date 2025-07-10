#!/bin/bash


ref_fa="/public/home/shengli-lab/shengnanluo/reference/star/human/v41/gencode.v41.GRCh38.primary_assembly.genome.fa" 
input_dir="/public/nanopore_cDNA/fastq"
output_dir="/public/nanopore_cDNA/bam"
mkdir -p ${output_dir}


samples=($(ls ${input_dir}))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}  


in_fq="${input_dir}/${sample}/filtered_${sample}_nanofilt.fastq"
out_sam="${output_dir}/${sample}.sam"
out_bam_tmp="${output_dir}/${sample}.tmp.bam"
out_bam="${output_dir}/${sample}.sorted.bam"


minimap2 -ax splice ${ref_fa} ${in_fq} -t 6 -o ${out_sam} && \
samtools view -bS ${out_sam} --threads 6 > ${out_bam_tmp} && \
samtools sort ${out_bam_tmp} --threads 6 -o ${out_bam} && \
samtools index ${out_bam}