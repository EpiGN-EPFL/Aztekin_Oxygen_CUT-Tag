import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def read_total_counts(indir, frip_file):
    """Reads the FRiP total counts file."""
    df = pd.read_csv(os.path.join(indir, frip_file), sep='\t', index_col=0)
    df.columns = [os.path.basename(col).removesuffix('.bam') for col in df.columns]
    return df

def read_spikein_counts(indir, spikein_file, total_columns):
    """Reads the spike-in counts file and aligns it with total count columns."""
    df_spikein = pd.read_csv(os.path.join(indir, spikein_file), sep='\t', quotechar="'")
    df_spikein.columns = [os.path.basename(col).removesuffix('.filtered.bam') for col in df_spikein.columns]
    count_spikein = df_spikein.sum(axis=0, numeric_only=True)
    count_spikein = count_spikein.reindex(total_columns, fill_value=0)
    return count_spikein

def calculate_spikein_ratios(df_total, count_spikein):
    """Appends spike-in counts and calculates percentage metrics."""
    df = pd.concat([df_total, pd.DataFrame([count_spikein], index=["SpikeIn"])])
    df.loc['SpikeIn_percentOverTotal'] = df.loc['SpikeIn'] / df.loc['Total_Read']
    df.loc['SpikeIn_percentOverAssigned'] = df.loc['SpikeIn'] / df.loc['Assigned']
    return df

def custom_sort_key(filename, sort_priority_cond, sort_priority_histone):
    """Defines custom sorting logic for filenames."""
    histone_priority = len(sort_priority_histone)
    condition_priority = len(sort_priority_cond)
    
    for histone_key in sort_priority_histone:
        if histone_key in filename:
            histone_priority = sort_priority_histone[histone_key]
            break
    for cond_key in sort_priority_cond:
        if cond_key in filename:
            condition_priority = sort_priority_cond[cond_key]
            break
    return (histone_priority, condition_priority, filename)

def sort_columns(df, sort_priority_cond, sort_priority_histone):
    """Sorts DataFrame columns using custom sort logic."""
    return df[sorted(df.columns, key=lambda col: custom_sort_key(col, sort_priority_cond, sort_priority_histone))]

def plot_data_and_save(df, cols, name, indir, output_file):
    """Plots data and saves the figure."""
    plot_data = df.loc[cols].transpose().reset_index()
    plot_data = plot_data.melt(id_vars='index', var_name='Status', value_name=name)

    plt.figure(figsize=(12, 6))
    sns.lineplot(data=plot_data, x='index', y=name, hue='Status', marker='o')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(indir, output_file), bbox_inches='tight')

def main(args):
    # Read input files
    df_total = read_total_counts(args.indir, args.frip)
    df_total = df_total.loc[:, ~df_total.columns.str.contains(r'(H3K27ac|H3K9me3)')]
    count_spikein = read_spikein_counts(args.indir, args.spikein, df_total.columns)

    # Process and sort data
    df = calculate_spikein_ratios(df_total, count_spikein)
    df = sort_columns(df, args.sort_priority_cond, args.sort_priority_histone)

    # Define plot parameters
    plots = [
        {
            'cols': ['Assigned', 'SpikeIn', 'Total_Read'],
            'name': 'Count',
            'output_file': 'Count.pdf'
        },
        {
            'cols': ['SpikeIn_percentOverTotal', 'SpikeIn_percentOverAssigned'],
            'name': 'Percentage',
            'output_file': 'SpikeIn_Percentage.pdf'
        },
        {
            'cols': ['FRiP'],  
            'name': 'Percentage',
            'output_file': 'FRiP.pdf'
        }
    ]

    # Loop through plot configurations and generate plots
    for plot_config in plots:
        plot_data_and_save(
            df=df,
            cols=plot_config['cols'],
            name=plot_config['name'],
            indir=args.indir,
            output_file=plot_config['output_file']
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spike-in and total counts data.")
    parser.add_argument('-i', '--indir', type=str, required=True)
    parser.add_argument('--spikein', type=str, required=False, default='multiBamSummary.spike_in.rawcounts.tsv', help="Path to spike-in counts file.")
    parser.add_argument('--frip', type=str, required=False, default='FRiP.csv', help="Path to FRiP counts file.")
    # parser.add_argument('--output', type=str, required=True, help="Path to save the output plot.")
    parser.add_argument('--sort_priority_cond', type=dict, default={"ALI": 0, "96-well": 1}, help="Sorting priority for conditions.")
    parser.add_argument('--sort_priority_histone', type=dict, default={
        "H3K27me3": 0, "H3K27ac": 1, "H3K4me3": 2, "H3K9me3": 3, "H3_": 4
    }, help="Sorting priority for histone modifications.")

    args = parser.parse_args()
    main(args)