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

input_dir = "/public/nanopore_cDNA/fastq"
output_csv = "/public/nanopore_cDNA/data/nanopore_read_stats.csv"

def process_sample(sample_dir):
    #Processes a single sample directory to extract read length and average quality
    dir_path = os.path.join(input_dir, sample_dir)
    
    target_file = f"filtered_{sample_dir}_nanofilt.fastq"
    fastq_path = os.path.join(dir_path, target_file)
    
    if not os.path.isfile(fastq_path):
        print(f"Warning: {target_file} not found in {dir_path}")
        return []
    
    print(f"Processing: {fastq_path}")
    results = []

    # Parse the FASTQ filec
    with open(fastq_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # Parse the FASTQ file
            length = len(record.seq)

            # Calculate average quality scorec
            qualities = record.letter_annotations["phred_quality"]
            avg_q = sum(qualities) / len(qualities)
            
            results.append([sample_dir, length, avg_q])
    
    return results

if __name__ == "__main__":
    # Get a list of all sample directories in the input directory
    sample_dirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]

    # Set the number of processes for the poolc
    num_processes = min(cpu_count(), 20)

    # Set the number of processes for the pool
    with Pool(processes=num_processes) as pool:
        all_results = pool.map(process_sample, sample_dirs)

    # Write all results to a single CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header rowc
        writer.writerow(['sample', 'read_length', 'avg_quality'])

        # Write header row
        for results in all_results:
            if results:  
                writer.writerows(results)
    
    print(f"Data processing completed. Results saved to: {output_csv}")
