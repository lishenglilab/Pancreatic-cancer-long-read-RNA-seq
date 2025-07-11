import pysam
import glob
from multiprocessing import Pool, cpu_count

def deduplicate_bam(bam_file):
    output_file = bam_file.replace(".duplicates_marked.bam", ".deduplicated.bam")
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam_in, \
             pysam.AlignmentFile(output_file, "wb", template=bam_in) as bam_out:
            
            for read in bam_in:
                if not read.is_duplicate:
                    bam_out.write(read)
        
        print(f"Completed: {bam_file} -> {output_file}")
        return True
    
    except Exception as e:
        print(f"Failed to process: {bam_file} | Error: {str(e)}")
        return False

if __name__ == "__main__":
    bam_files = glob.glob("*.duplicates_marked.bam")
    print(f"Found {len(bam_files)} files to process")
    
    num_processes = min(20, cpu_count())
    with Pool(processes=num_processes) as pool:
        results = pool.map(deduplicate_bam, bam_files)
    
    success_count = sum(results)
    print(f"Processing completed: {success_count}/{len(bam_files)} succeeded")
