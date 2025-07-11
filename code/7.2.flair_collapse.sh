#!/bin/bash


ref_gtf="/public/home/shengli-lab/shengnanluo/reference/star/human/v41/gencode.v41.primary_assembly.annotation.gtf"
ref_fa="/public/home/shengli-lab/shengnanluo/reference/star/human/v41/gencode.v41.GRCh38.primary_assembly.genome.fa"
bed_file="/public/nanopore_cDNA/filar/ffilar/merged_corrected.bed"
output_dir="/public/nanopore_cDNA/filar/ffilar"
fastq_dir="/public/nanopore_cDNA/fastq"

mkdir -p "${output_dir}"
mkdir -p "${output_dir}/tmp3"


fastq_files=()
for sample in $(ls "${fastq_dir}"); do
    fastq_files+=("${fastq_dir}/${sample}/filtered_${sample}_nanofilt.fastq")
done

# 合并FASTQ路径(更安全的方式)
fastq_list=$(printf "%s," "${fastq_files[@]}" | sed 's/,$//')

# 输出合并后的FASTQ列表（调试用）
echo "所有FASTQ文件: ${fastq_list}"

# 运行flair collapse
flair collapse -g "${ref_fa}" -q "${bed_file}" --gtf "${ref_gtf}" \
     -r "${fastq_list}" \
     --output "${output_dir}/all_samples_collapsed3" \
     --threads 50 \
     --temp_dir "${output_dir}/tmp3"

echo "所有样本处理完成"
