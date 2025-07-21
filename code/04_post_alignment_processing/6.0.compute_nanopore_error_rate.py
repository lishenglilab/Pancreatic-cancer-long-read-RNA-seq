"""
This script calculates the overall alignment error rate for a nanopore sequencing BAM file.
It parses the CIGAR string and the NM (edit distance) tag for each read to quantify
matches, mismatches, insertions, and deletions, providing a summary of alignment accuracy.
"""
import pysam

def compute_error_rate(bam_path):
    """
    Computes and prints the error rate statistics for a given BAM file.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Initialize counters for overall statistics
    total_matches = 0
    total_mismatches = 0
    total_insertions = 0
    total_deletions = 0

    for read in bam.fetch(until_eof=True):
        # Skip alignments that are not primary or are unmapped
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Initialize statistics for the current read
        read_matches = 0
        read_insertions = 0
        read_deletions = 0

        # Tally operations from the CIGAR string
        for op, length in read.cigartuples:
            if op == 0:  
                read_matches += length
            elif op == 1: 
                read_insertions += length
            elif op == 2:  
                read_deletions += length

        # The 'NM' tag represents the total edit distance (mismatches + insertions + deletions).
        try:
            nm = read.get_tag('NM')
            # NM = mismatches + insertions + deletions
            read_mismatches = nm - read_insertions - read_deletions
            if read_mismatches < 0:
                print(f"Warning: Negative mismatches for read {read.query_name}. NM={nm}, I={read_insertions}, D={read_deletions}")
                # Removed the line that sets read_mismatches = 0
        except KeyError:
            read_mismatches = 0  

        # Add the current read's stats to the overall total
        total_matches += read_matches
        total_mismatches += read_mismatches
        total_insertions += read_insertions
        total_deletions += read_deletions

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
