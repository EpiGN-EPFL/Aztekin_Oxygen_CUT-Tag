import numpy as np
import pandas as pd
import os
import glob
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import re

def replace_dot(value):
    """Convert '.' to NaN, otherwise convert to float."""
    return np.nan if value == '.' else float(value)

colors = {
    "H3K27me3": "#027BB8",
    "H3K27ac": "#009697",
    "H3K4me3": "#E56F8E",
    "H3K9me3":"#a63c8d",
    "H3": "#808080"
}

# Sort df columns 
sort_priority_histone = {
    "H3K27me3": 0,
    "H3K27ac": 1,
    "H3K4me3":2,
    "H3K9me3":3,
    "H3_":4
}

def custom_sort_key(filename):
    # Set default large values to move unmatched items to the end
    histone_priority = len(sort_priority_histone)
    
    # Determine the histone modification sort order
    for histone_key in sort_priority_histone:
        if histone_key in filename:
            histone_priority = sort_priority_histone[histone_key]
            break

    # Return a tuple that sorts by histone modification first, then filename
    return (histone_priority, filename)

# Define a function to extract components from column names
def extract_components(col_name):
    pattern = r'^(.*?)_(.*?)_(.*?)_(.*?)_(.*?)_(.*?)$'  # Regex for the pattern
    match = re.match(pattern, col_name)
    if match:
        return [
            match.group(1),  # date
            match.group(2),  # initial
            match.group(3),  # condition
            match.group(4),  # histone
            match.group(5),  # operation
            match.group(6),  # calculation
        ]
    return [None] * 6

def read_file(file):
    df = pd.read_csv(file, sep='\t', nrows=2)
    cols = [col for col in df.columns if 'mean' in col and 'H3K27ac' not in col]
    df =pd.read_csv(file, sep='\t', usecols=cols + ['chr'], converters={col: replace_dot for col in cols})
    #df=df[~df['chr'].str.contains('S')]
    df=df.drop(columns=['chr'])
    return (df)

def parse_args():
    # Create the parser
    parser = argparse.ArgumentParser(description="Violin plot for bed files (tsv format) in a folder.")

    # Add the argument
    parser.add_argument('-i', '--input_dir', type=str, required=True, 
                        help="The path to the directory containing files to read.")
    parser.add_argument('-p', '--pattern', type=str, required=False, 
                        default="_combined.bed", help="Patten of the files to read in input")
    parser.add_argument('-v', '--verbose', type=bool, required=False, default=False)
    
    # Parse the arguments
    args = parser.parse_args()

    return(args)

def main():
    # Read summary from featureCounts
    args=parse_args()
    files = glob.glob(os.path.join(args.input_dir, f'*{args.pattern}'))

    for file in files:
        print(f'reading file {file}')
        df = read_file(file)
        
        # Calculate the median of each column
        column_medians = df.median()  # Adjust the row slicing (6:) to exclude metadata row

        # Extract the components 
        components = [extract_components(col) for col in df.columns]
        component_labels = ['date', 'initial', 'condition', 'histone', 'operation', 'calculation']
        components_df = pd.DataFrame(components, columns=component_labels, index=df.columns)

        # Combine components and medians into a single DataFrame for plotting
        components_df['median'] = column_medians.values  # Align medians with components_df

        # Alternatively, plot a box plot
        plt.figure(figsize=(8, 6))
        sns.boxplot(data=components_df, x='histone', y='median', hue='histone',palette=colors)
        sns.stripplot(data=components_df, x='histone', y='median', alpha=0.6, hue='date')  # Overlay dots for clarity
        plt.ylim(-2.5, 2.5)
        plt.axhline(0, color='gray', linestyle='--', linewidth=1)
        plt.title(f'{os.path.basename(file)}, n(L+S)={df.shape[0]}')
        plt.xlabel('')
        plt.ylabel(", ".join(components_df['condition'].unique()))
        savepath=f'{file.removesuffix(args.pattern)}.boxplot.pdf'
        plt.savefig(savepath, bbox_inches='tight')
        plt.close()
        print(f"Saving result to {savepath}")

if __name__ == "__main__":
    main()