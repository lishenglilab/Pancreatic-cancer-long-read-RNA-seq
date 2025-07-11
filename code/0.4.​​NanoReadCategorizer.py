import os
import csv
from multiprocessing import Pool

def parse_fastq(file_path):
    with open(file_path, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header: break
            seq = f.readline().strip()
            _ = f.readline() 
            qual = f.readline().strip()
            if len(seq) != len(qual):
                print(f"Format error: sequence and quality length mismatch in {file_path}")
                continue
            yield (seq, qual)

def calculate_metrics(seq, qual):
    read_length = len(seq)
    qual_scores = [ord(c)-33 for c in qual]
    avg_qual = sum(qual_scores)/len(qual_scores) if qual_scores else 0
    return read_length, avg_qual

def process_single_file(file_path):
    counts = [0]*4 
    try:
        for seq, qual in parse_fastq(file_path):
            length, qual = calculate_metrics(seq, qual)
            
            if length < 200:
                if qual < 7: counts[0] += 1
                else: counts[1] += 1
            else:
                if qual < 7: counts[2] += 1
                else: counts[3] += 1
        
        total = sum(counts)
        ratios = [f"{(c/total*100):.2f}%" if total else "0.00%" for c in counts]
        return {
            "filename": os.path.basename(file_path),
            "condition1_count": counts[0],
            "condition1_ratio": ratios[0],
            "condition2_count": counts[1],
            "condition2_ratio": ratios[1],
            "condition3_count": counts[2],
            "condition3_ratio": ratios[2],
            "condition4_count": counts[3],
            "condition4_ratio": ratios[3],
            "total_reads": total
        }
    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}")
        return None

def main():
    input_root = "/public/nanopore_cDNA/fastq"
    output_dir = "/public/nanopore_cDNA/data/"
    os.makedirs(output_dir, exist_ok=True)
    
    file_list = []
    for dir_name in os.listdir(input_root):
        dir_path = os.path.join(input_root, dir_name)
        if os.path.isdir(dir_path):
            fastq_file = os.path.join(dir_path, f"filtered_{dir_name}.fastq")
            if os.path.exists(fastq_file):
                file_list.append(fastq_file)
            else:
                print(f"Warning: {fastq_file} not found")

    with Pool(20) as pool:
        results = [res for res in pool.map(process_single_file, file_list) if res]

    output_path = os.path.join(output_dir, "nanopore_qc_report.csv")
    fieldnames = [
        "filename", "total_reads",
        "condition1_count", "condition1_ratio",
        "condition2_count", "condition2_ratio",
        "condition3_count", "condition3_ratio",
        "condition4_count", "condition4_ratio"
    ]
    
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Analysis completed. Results saved to {output_path}")

if __name__ == "__main__":
    main()
