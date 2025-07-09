import pysam

def compute_error_rate(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")

    total_matches = 0
    total_mismatches = 0
    total_insertions = 0
    total_deletions = 0

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # 初始化当前read的统计
        read_matches = 0
        read_insertions = 0
        read_deletions = 0

        # 统计CIGAR操作（当前read）
        for op, length in read.cigartuples:
            if op == 0:  # M（匹配或错配）
                read_matches += length
            elif op == 1:  # I（插入）
                read_insertions += length
            elif op == 2:  # D（缺失）
                read_deletions += length

        # 通过NM标签计算当前read的错配数
        try:
            nm = read.get_tag('NM')
            # NM = mismatches + insertions + deletions
            read_mismatches = nm - read_insertions - read_deletions
            if read_mismatches < 0:
                print(f"Warning: Negative mismatches for read {read.query_name}. NM={nm}, I={read_insertions}, D={read_deletions}")
                # Removed the line that sets read_mismatches = 0
        except KeyError:
            read_mismatches = 0  # 无NM标签则跳过

        # 累加到全局统计
        total_matches += read_matches
        total_mismatches += read_mismatches
        total_insertions += read_insertions
        total_deletions += read_deletions

    # 计算总错误率和输出
    total_errors = total_mismatches + total_insertions + total_deletions
    total_aligned = total_matches + total_errors
    error_rate = total_errors / total_aligned if total_aligned > 0 else 0

    print(f"Total Matches     : {total_matches}")
    print(f"Total Mismatches  : {total_mismatches}")
    print(f"Total Insertions  : {total_insertions}")
    print(f"Total Deletions   : {total_deletions}")
    print(f"Total Aligned     : {total_aligned}")
    print(f"Total Errors      : {total_errors}")
    print(f"Error Rate        : {error_rate:.4f} ({error_rate*100:.2f}%)")

    bam.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python compute_nanopore_error_rate.py aligned.bam")
        sys.exit(1)
    compute_error_rate(sys.argv[1])
