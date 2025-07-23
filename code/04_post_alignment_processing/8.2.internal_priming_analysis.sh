import os
import multiprocessing as mp
from collections import defaultdict
import pandas as pd

def parse_fasta(filename):
    """Parses FASTA files and returns a list of DNA sequences (uppercase)"""
    sequences = []
    current_seq = ""

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq.upper())
                    current_seq = ""
            else:
                current_seq += line

        if current_seq:
            sequences.append(current_seq.upper())

    return sequences

def analyze_single_file(filename):
    """Analyzes a single FASTA file and calculates key statistics"""
    try:
        sequences = parse_fasta(filename)
        total_sequences = len(sequences)

        # Handle empty files
        if total_sequences == 0:
            return {
                'filename': filename,
                'total_sequences': 0,
                'avg_a_percentage': 0,
                'sequences_with_6plus_a': 0,
                'percentage_with_6plus_a': 0,
                'individual_a_percentages': []
            }

        a_percentages = []
        sequences_with_6plus_a = 0

        for seq in sequences:
            if seq:  # Only process non-empty sequences
                # Calculate A content percentage
                a_count = seq.count('A')
                a_percentage = (a_count / len(seq)) * 100
                a_percentages.append(a_percentage)
                
                # Check for ≥6 consecutive As (poly-A stretch)
                if 'AAAAAA' in seq:
                    sequences_with_6plus_a += 1

        # Calculate summary statistics
        avg_a_percentage = sum(a_percentages) / len(a_percentages)
        percentage_with_6plus_a = (sequences_with_6plus_a / total_sequences) * 100

        return {
            'filename': filename,
            'total_sequences': total_sequences,
            'avg_a_percentage': avg_a_percentage,
            'sequences_with_6plus_a': sequences_with_6plus_a,
            'percentage_with_6plus_a': percentage_with_6plus_a,
            'individual_a_percentages': a_percentages
        }

    except Exception as e:
        print(f"Error processing {filename}: {e}")
        return None

def main():
    # List of FASTA files to process (paired samples)
    fasta_files = [
        'APSC-1-1_flanks.fa', 'APSC-1-2_flanks.fa',
        'BXPC3-1_flanks.fa', 'BXPC3-2_flanks.fa',
        'Capan-1-1_flanks.fa', 'Capan-1-2_flanks.fa',
        'Capan-2-1_flanks.fa', 'Capan-2-2_flanks.fa',
        'HuP-T4-1_flanks.fa', 'HuP-T4-2_flanks.fa',
        'Mia-Paca-2-1_flanks.fa', 'Mia-Paca-2-2_flanks.fa',
        'PANC0203-1_flanks.fa', 'PANC0203-2_flanks.fa',
        'Panc0327-1_flanks.fa', 'Panc0327-2_flanks.fa',
        'PANC1005-1_flanks.fa', 'PANC1005-2_flanks.fa',
        'SW1990-1_flanks.fa', 'SW1990-2_flanks.fa'
    ]

    # Process files in parallel (20 processes)
    with mp.Pool(processes=20) as pool:
        results = pool.map(analyze_single_file, fasta_files)

    # Filter out failed analyses
    results = [r for r in results if r is not None]

    # Generate summary report
    print("="*80)
    print("FASTA File A-Base Analysis Report")
    print("="*80)

    # Create results table
    summary_data = []
    for result in results:
        summary_data.append({
            'Filename': result['filename'],
            'Total Sequences': result['total_sequences'],
            'Avg A %': f"{result['avg_a_percentage']:.2f}",
            '≥6A Sequences': result['sequences_with_6plus_a'],
            '% with ≥6A': f"{result['percentage_with_6plus_a']:.2f}"
        })

    # Format and display table
    df = pd.DataFrame(summary_data)
    print(df.to_string(index=False))

    # Calculate overall statistics
    total_sequences = sum(r['total_sequences'] for r in results)
    total_with_6plus_a = sum(r['sequences_with_6plus_a'] for r in results)
    overall_6plus_a_percentage = (total_with_6plus_a / total_sequences) * 100
    
    # Calculate global A percentage average
    all_a_percentages = []
    for result in results:
        all_a_percentages.extend(result['individual_a_percentages'])
    overall_avg_a = sum(all_a_percentages) / len(all_a_percentages)

    # Print aggregate stats
    print("\n" + "="*80)
    print("Aggregate Statistics")
    print("="*80)
    print(f"Total Sequences Analyzed: {total_sequences:,}")
    print(f"Global Average A Percentage: {overall_avg_a:.2f}%")
    print(f"Sequences with ≥6 Consecutive As: {total_with_6plus_a:,}")
    print(f"Percentage with Poly-A Stretches: {overall_6plus_a_percentage:.2f}%")

    # Save detailed results
    df.to_csv('fasta_analysis_results.csv', index=False, encoding='utf-8-sig')
    print(f"\nResults saved to: fasta_analysis_results.csv")

if __name__ == "__main__":
    main()