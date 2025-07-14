"""
This script generates a read length distribution plot for multiple samples.
It calculates the kernel density estimation (KDE) of read lengths for each sample
and plots them on a single graph. It also marks the median peak position across all samples.
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from multiprocessing import Pool

folders = [
    '02.03-1', '02.03-2', '1005-1', '1005-2',
    '3.27-1', '3.27-2', 'APSC-1-1', 'APSC-1-2',
    'BXPC3-1', 'BXPC3-2', 'Capan-1-1', 'Capan-1-2',
    'Capan2-1', 'Capan2-2', 'HupT4-1', 'HupT4-2',
    'Mia-Paca-1', 'Mia-Paca-2', 'SW1990-1', 'SW1990-2'
]


def load_data(folder):
    #Load read length data from a CSV file for a given sample folder.
    try:
        file_path = os.path.join(folder, 'csv_data', 'length.csv')
        df = pd.read_csv(file_path)
        data = df['ReadLength'].dropna().values
        return (folder, data)
    except Exception as e:
        print(f"Error processing {folder}: {str(e)}")
        return (folder, None)


if __name__ == '__main__':
    # Use multiprocessing to load data in parallel
    with Pool(processes=20) as pool:
        results = pool.map(load_data, folders)

    # Process results from data loading
    data_dict = {}
    all_data = []
    for folder, data in results:
        if data is not None and len(data) > 0:
            data_dict[folder] = data
            all_data.extend(data)

    all_data = np.array(all_data)

    print(f"[Debug] Original data range: {np.min(all_data):.1f} - {np.max(all_data):.1f} (total {len(all_data)} data points)")

    # Define percentile range to filter out outliers
    lower_percentile = 1  
    upper_percentile = 99  

    # Calculate percentile values
    q_low = np.percentile(all_data, lower_percentile)
    q_high = np.percentile(all_data, upper_percentile)

    print(f"[Debug] Using quantile range: {q_low:.1f} - {q_high:.1f} ({lower_percentile}%–{upper_percentile}% quantiles)")
    print(f"[Debug] Left tail removed (%): {np.mean(all_data < q_low) * 100:.2f}%")
    print(f"[Debug] Right tail removed (%): {np.mean(all_data > q_high) * 100:.2f}%")

    # Generate x-axis values for the plot
    x = np.linspace(q_low, q_high, 1000)

    # Set up the plot
    plt.figure(figsize=(15, 8))
    plt.rcParams.update({
        'font.size': 12,
        'pdf.fonttype': 42,
        'axes.axisbelow': True 
    })
    colors = plt.cm.tab20.colors

    peak_x_values = []


    def compute_kde(args):
        #Compute Kernel Density Estimation for a given sample's data
        folder, data = args
        if len(data) < 2:
            return None
        kde = gaussian_kde(data)
        y = kde(x)
        y_count = y * len(data)
        return (folder, x, y_count)


    # Compute KDE for each sample in parallel
    with Pool(processes=20) as pool:
        kde_results = pool.map(compute_kde, [(k, v) for k, v in data_dict.items()])

    # Plot KDE results for each sample
    for i, result in enumerate(kde_results):
        if result is not None:
            folder, x_vals, y_vals = result
            plt.plot(x_vals, y_vals, color=colors[i], label=folder, alpha=0.8, linewidth=1.5)

            # Find and store the peak of the KDE curve
            peak_index = np.argmax(y_vals)
            peak_x = x_vals[peak_index]
            peak_x_values.append(peak_x)

    # Calculate and plot the median peak position
    if peak_x_values:
        median_peak = np.median(peak_x_values)
        print(f"[Debug] Median peak position across samples: {median_peak:.1f} bp")

        plt.axvline(median_peak,
                    color='black',
                    linestyle='--',
                    linewidth=2,
                    alpha=0.9,
                    zorder=100)

        plt.text(median_peak, plt.ylim()[1] * 0.85,
                 f'Median Peak Position\n{median_peak:.0f} bp',
                 ha='center',
                 va='top',
                 fontsize=10,
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    # Set plot labels and title
    plt.xlabel('Read Length (bp)', fontsize=14)
    plt.ylabel('Density-Scaled Count', fontsize=14)
    plt.title(f'Read Length Distribution ({lower_percentile}%-{upper_percentile}% Quantile Range)', fontsize=16, pad=20)

    # Configure the legend
    legend = plt.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.,
                        frameon=False,
                        fontsize=10,
                        title='Sample Names',
                        title_fontsize=12)

    # Customize plot aesthetics
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle='--', which='both')

    # Add lines for quantile boundaries
    ax.axvline(q_low, color='gray', linestyle=':', alpha=0.6, linewidth=1)
    ax.axvline(q_high, color='gray', linestyle=':', alpha=0.6, linewidth=1)

    # Add text for quantile boundaries
    plt.text(q_low, ax.get_ylim()[1] * 0.95,
             f'{lower_percentile}%: {q_low:.0f} bp',
             ha='left', va='top', alpha=0.7)
    plt.text(q_high, ax.get_ylim()[1] * 0.95,
             f'{upper_percentile}%: {q_high:.0f} bp',
             ha='right', va='top', alpha=0.7)

    # Save the plot to a PDF file
    plt.savefig('Read_Length_Distribution.pdf',
                format='pdf',
                dpi=600,
                bbox_inches='tight',
                facecolor='white',
                edgecolor='none',
                transparent=False)

    plt.close()
    print("[Status] Chart saved as Read_Length_Distribution.pdf")
