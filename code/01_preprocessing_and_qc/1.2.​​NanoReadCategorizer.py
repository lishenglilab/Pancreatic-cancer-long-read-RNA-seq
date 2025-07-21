"""
This script processes filtered FASTQ files to categorize reads based on length and quality.
It classifies each read into one of four conditions
summarizing the counts and ratios for each category per sample file.
The conditions are:
1. Length < 200 and Quality < 7
2. Length < 200 and Quality >= 7
3. Length >= 200 and Quality < 7
4. Length >= 200 and Quality >= 7
"""
import pandas as pd
import multiprocessing as mp
import numpy as np
from functools import reduce
import os

def classify_and_count(df):
    """
    Classifies a chunk of a DataFrame based on length and quality,
    and returns the counts for each sample.
    """
    # Dictionary to store the results for this chunk
    results = {}
    # Group by the 'sample' column
    for sample, group in df.groupby('sample'):
        # Initialize counters for the current sample
        counts = {
            'Condition 1 (length < 200, quality < 7)': 0,
            'Condition 2 (length < 200, quality >= 7)': 0,
            'Condition 3 (length >= 200, quality < 7)': 0,
            'Condition 4 (length >= 200, quality >= 7)': 0
        }
        # Iterate over the rows in the group and classify them
        for _, row in group.iterrows():
            length = row['read_length']
            quality = row['avg_quality_correct']
            if length < 200 and quality < 7:
                counts['Condition 1 (length < 200, quality < 7)'] += 1
            elif length < 200 and quality >= 7:
                counts['Condition 2 (length < 200, quality >= 7)'] += 1
            elif length >= 200 and quality < 7:
                counts['Condition 3 (length >= 200, quality < 7)'] += 1
            elif length >= 200 and quality >= 7:
                counts['Condition 4 (length >= 200, quality >= 7)'] += 1
        results[sample] = counts
    return results

def merge_results(dict1, dict2):
    """
    Merges two count dictionaries, summing the values for common keys.
    """
    for sample, counts in dict2.items():
        if sample in dict1:
            for condition, count in counts.items():
                dict1[sample][condition] += count
        else:
            dict1[sample] = counts
    return dict1

def main():
    # CSV file path - update with your filename
    file_path = 'nanopore_read_stats_correct_avg_q_AfterKraken.csv'

    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found. Please ensure the file path is correct.")
        return

    # Number of CPU cores to use
    num_processes = 40

    print(f"Reading data from '{file_path}'...")
    # Read the data
    df = pd.read_csv(file_path)
    print("Data reading complete.")

    # Split the DataFrame into chunks for parallel processing
    df_split = np.array_split(df, num_processes)

    print(f"Starting parallel computation with {num_processes} CPU cores...")
    # Create a pool of worker processes
    with mp.Pool(num_processes) as pool:
        # Process all data chunks in parallel
        results = pool.map(classify_and_count, df_split)

    print("Computation finished, merging results...")
    # Merge the results returned by all processes
    # We need to filter out any dictionaries that might be empty
    final_counts = reduce(merge_results, [res for res in results if res])
    print("Results merged.")

    # --- New code: Calculate the proportion of reads for each category in each sample ---
    final_proportions = {}
    for sample, counts in final_counts.items():
        total_reads = sum(counts.values())
        proportions = {}
        if total_reads > 0:
            for condition, count in counts.items():
                proportions[condition] = (count / total_reads) * 100
        else:
            for condition, count in counts.items():
                proportions[condition] = 0
        final_proportions[sample] = proportions
    # --- End of new code ---

    # --- Modified code: Print the final proportion results ---
    print("\n--- Classification and Proportion Results ---")
    for sample, proportions in sorted(final_proportions.items()):
        print(f"\nSample: {sample}")
        # Get and print the total number of reads for each sample
        total_reads = sum(final_counts[sample].values())
        print(f"  Total reads: {total_reads}")
        for condition, proportion in proportions.items():
            # Get the original count to display it as well
            count = final_counts[sample][condition]
            print(f"  {condition}: {proportion:.2f}% ({count} reads)")
    print("\n----------------------")
    # --- End of modified code ---

if __name__ == '__main__':
    main()
if __name__ == "__main__":
    main()
