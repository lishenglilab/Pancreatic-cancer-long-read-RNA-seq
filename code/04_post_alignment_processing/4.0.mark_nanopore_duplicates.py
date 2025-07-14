"""
This script marks PCR duplicates in nanopore sequencing data (BAM format). This script
adopts a more robust strategy by first clustering reads based on approximate genomic
location and splice junctions, and then comparing sequence similarity within each
cluster to identify and mark duplicates.
"""

import pysam
import Levenshtein
from collections import defaultdict

# Parameters (can adjust these as needed)
POSITION_TOLERANCE = 20
SEQUENCE_SIMILARITY_THRESHOLD = 0.95

def get_junction_signature(read):
    """
    Extracts a signature of splice junctions from a read's CIGAR string.c
    """
    junctions = []
    pos = read.reference_start
    for cigartuple in read.cigartuples:
        op, length = cigartuple
        if op == 3:  # splice junction ("N" in CIGAR)
            junctions.append((pos, pos + length))
        if op in (0, 2, 3):  # M, D, N
            pos += length
    return tuple(junctions)

def cluster_reads(bamfile):
    """
    Clusters reads from a BAM file based on genomic location and splice junctions.
    """
    clusters = defaultdict(list)

    for read in bamfile.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        junctions = get_junction_signature(read)
        cluster_key = (
            read.reference_name,
            round(read.reference_start / POSITION_TOLERANCE),
            round(read.reference_end / POSITION_TOLERANCE),
            junctions
        )
        clusters[cluster_key].append(read)

    return clusters

def sequence_similarity(seq1, seq2):
    """
    Uses Levenshtein distance to measure sequence differences
    """
    edit_distance = Levenshtein.distance(seq1, seq2)
    max_len = max(len(seq1), len(seq2))
    similarity = 1 - (edit_distance / max_len)
    return similarity

def mark_duplicates(input_bam, output_bam):
    """
    - Opens input/output BAM files
    - Clusters reads by genomic location and splice junctions
    - Within each cluster: keeps first read as representative
    - Marks subsequent reads as duplicates if sequence similarity > threshold
    - Writes all reads to output (marked duplicates still included)
    """
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    clusters = cluster_reads(bam_in)
    duplicates_marked = 0

    for reads in clusters.values():
        if len(reads) == 1:
            bam_out.write(reads[0])
            continue

        representative = reads[0]
        bam_out.write(representative)

        for read in reads[1:]:
            sim = sequence_similarity(representative.query_sequence, read.query_sequence)
            if sim >= SEQUENCE_SIMILARITY_THRESHOLD:
                read.is_duplicate = True
                duplicates_marked += 1
            bam_out.write(read)

    bam_in.close()
    bam_out.close()

    print(f"Duplicates marked: {duplicates_marked}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python mark_nanopore_duplicates.py input.bam output_marked.bam")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    mark_duplicates(input_bam, output_bam)
