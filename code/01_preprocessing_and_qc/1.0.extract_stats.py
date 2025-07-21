"""
This script processes FASTQ files from multiple sample directories to extract read statistics.
For each read in a filtered FASTQ file, it calculates the read length and the average quality score.
The results (sample, read_length, avg_quality) are then compiled and saved into a single CSV file.
It uses multiprocessing to process samples in parallel for efficiency.
"""
import os
from Bio import SeqIO
import csv
from multiprocessing import Pool, cpu_count
import math  # Import math library for log calculations

# Set input and output paths
input_dir = "/public/nanopore_cDNA/fastq"
output_csv = "/public/nanopore_cDNA/data/nanopore_read_stats_correct_avg_q.csv"

# Define function to process a single sample (using correct Q-value averaging method)
def process_sample(sample_dir):
    dir_path = os.path.join(input_dir, sample_dir)
    target_file = f"filtered_{sample_dir}_nanofilt.fastq"
    fastq_path = os.path.join(dir_path, target_file)

    if not os.path.isfile(fastq_path):
        print(f"Warning: {target_file} not found in {dir_path}")
        return []

    print(f"Processing: {fastq_path}")
    results = []

    with open(fastq_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            length = len(record.seq)
            if length == 0:
                continue

            qualities = record.letter_annotations["phred_quality"]

            # --- Steps for correct average Q-value calculation ---
            # 1. Convert all Q-values to error probabilities P and sum them
            sum_of_error_probs = sum([10**(-q / 10.0) for q in qualities])

            # 2. Calculate average error probability P_avg
            avg_error_prob = sum_of_error_probs / length

            # 3. Convert average error probability back to final Q-value Q_avg
            # If avg_error_prob is 0 (meaning perfect quality), assign a high Q-value (e.g., 93)
            if avg_error_prob == 0:
                avg_q_correct = 93.0
            else:
                avg_q_correct = -10 * math.log10(avg_error_prob)
            # --- End of calculation ---

            results.append([sample_dir, length, avg_q_correct])

    return results

# Main program (same as before)
if __name__ == "__main__":
    sample_dirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    num_processes = min(cpu_count(), 20)

    with Pool(processes=num_processes) as pool:
        all_results = pool.map(process_sample, sample_dirs)

    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample', 'read_length', 'avg_quality_correct'])

        for results in all_results:
            if results:
                writer.writerows(results)

    print(f"Data processing complete. Results saved to: {output_csv}")
