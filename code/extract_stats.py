import os
from Bio import SeqIO
import csv
from multiprocessing import Pool, cpu_count

input_dir = "/public/nanopore_cDNA/fastq"
output_csv = "/public/nanopore_cDNA/data/nanopore_read_stats.csv"

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
            
            qualities = record.letter_annotations["phred_quality"]
            avg_q = sum(qualities) / len(qualities)
            
            results.append([sample_dir, length, avg_q])
    
    return results

if __name__ == "__main__":
    sample_dirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    
    num_processes = min(cpu_count(), 20)
    
    with Pool(processes=num_processes) as pool:
        all_results = pool.map(process_sample, sample_dirs)
    
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['sample', 'read_length', 'avg_quality'])
        
        for results in all_results:
            if results:  
                writer.writerows(results)
    
    print(f"Data processing completed. Results saved to: {output_csv}")
