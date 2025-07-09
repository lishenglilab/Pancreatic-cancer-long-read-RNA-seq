import os
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from matplotlib.ticker import FuncFormatter 



if __name__ == '__main__':
    folders = [
    '02.03-1', '02.03-2', '1005-1', '1005-2',
    '3.27-1', '3.27-2', 'APSC-1-1', 'APSC-1-2',
    'BXPC3-1', 'BXPC3-2', 'Capan-1-1', 'Capan-1-2',
    'Capan2-1', 'Capan2-2', 'HupT4-1', 'HupT4-2',
    'Mia-Paca-1', 'Mia-Paca-2', 'SW1990-1', 'SW1990-2'
]


def load_data(folder):
    try:
        file_path = os.path.join(folder, 'csv_data', 'quality.csv')
        df = pd.read_csv(file_path)
        data = df['AverageQuality'].dropna().values
        return (folder, data)
    except Exception as e:
        print(f"Error processing {folder}: {str(e)}")
        return (folder, None)


if __name__ == '__main__':
    with Pool(processes=20) as pool:
        results = pool.map(load_data, folders)

    data_dict = {}
    all_data = []
    for folder, data in results:
        if data is not None and len(data) > 0:
            data_dict[folder] = data
            all_data.extend(data)

    global_min = np.min(all_data)
    global_max = np.max(all_data)
    x = np.linspace(global_min, global_max, 1000)

    plt.figure(figsize=(15, 8), dpi=300)  
    colors = plt.cm.tab20.colors


    def compute_kde(args):
        folder, data = args
        if len(data) < 2:
            return None
        kde = gaussian_kde(data)
        y = kde(x)
        y_count = y * len(data)
        return (folder, x, y_count)


    with Pool(processes=20) as pool:
        kde_results = pool.map(compute_kde, [(k, v) for k, v in data_dict.items()])

    for i, result in enumerate(kde_results):
        if result is not None:
            folder, x_vals, y_vals = result
            plt.plot(x_vals, y_vals, color=colors[i],
                     label=folder, alpha=0.8, linewidth=1.5)


    plt.xlim(0, 40)  
    plt.ylim(0, 1.5e6)  

    plt.xlabel('Average Quality', fontsize=14)
    plt.ylabel('Count (×10⁶)', fontsize=14) 


    def million_formatter(x, pos):
        return f"{x / 1e6:.1f}"  


    plt.gca().yaxis.set_major_formatter(FuncFormatter(million_formatter))

    plt.margins(x=0.02, y=0.05)  

    plt.title('Quality Distribution (clean) ', fontsize=16, pad=20)
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.,
               fontsize=10)

    plt.tight_layout()
    output_path = 'quality_distribution_clean2.pdf' 
    plt.savefig(output_path, format='pdf', dpi=300,
                bbox_inches='tight', pad_inches=0.5)

    print(f"PDF saved to: {os.path.abspath(output_path)}")
    plt.close()