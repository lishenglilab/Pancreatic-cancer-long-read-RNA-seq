#!/bin/bash
#SBATCH --job-name=exon_coverage
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --account=shengli-lab

source /public/home/shengli-lab/shengnanluo/miniconda3/bin/activate pancreas

# ---------- 参数配置 ----------
GTF="/public/home/shengli-lab/shengnanluo/reference/star/human/v41/gencode.v41.primary_assembly.annotation.gtf"
BAM_DIR="/public/nanopore_cDNA/bam"
OUT_DIR="/public/nanopore_cDNA/bam/coverage_results"
mkdir -p $OUT_DIR

# ---------- 函数定义 ----------
process_bam() {
    local bam=$1
    local sample=$(basename "$bam" .dedup.sorted.bam)
    
    # 检查索引是否存在（兼容 .bam.bai 和 .bai 两种格式）
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        echo "[$(date)] Indexing $bam..."
        if ! samtools index "$bam"; then
            echo "ERROR: Failed to index $bam" >&2
            return 1
        fi
    fi

    # 计算覆盖度
    echo "[$(date)] Processing $sample..."
    bedtools coverage -a merged_exons.bed -b "$bam" -d > "$OUT_DIR/${sample}.coverage_per_base.txt" || {
        echo "ERROR: bedtools failed for $bam" >&2
        return 1
    }
}
export -f process_bam
export OUT_DIR

# ---------- 主流程 ----------
echo "=== 开始外显子覆盖度分析 ==="
date

# Step 1: 提取并合并外显子（仅运行一次）
if [[ ! -f "merged_exons.bed" ]]; then
    echo "[$(date)] Extracting exons from GTF..."
    grep -w "exon" "$GTF" | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5, ".", ".", $7}' > exons.bed || exit 1
    bedtools sort -i exons.bed | bedtools merge -s -i - > merged_exons.bed || exit 1
fi

# Step 2: 并行处理所有BAM文件
echo "[$(date)] Processing BAM files (20 parallel jobs)..."
find "$BAM_DIR" -name "*.dedup.sorted.bam" | parallel -j 10 --halt soon,fail=1 process_bam || {
    echo "ERROR: Parallel execution failed" >&2
    exit 1
}

# 完成
echo "[$(date)] All tasks completed successfully!"
echo "结果文件保存在: $OUT_DIR/"
