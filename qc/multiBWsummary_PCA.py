import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

file = '/scratch/hhu/aztekin_v2/analysis_v4/mouse/multiBigwigSummary_bamCoverage_scaleFactor_log2/PCA_multiBWsummary.10000bp.tsv'

# Read the file
df = pd.read_csv(file, sep='\t', skiprows=0, header=1)

# Drop the 'Eigenvalue' column 
df.drop(columns='Eigenvalue', inplace=True)

# Transpose the DataFrame so that the samples are the rows
df_t = df.set_index('Component').transpose()

pcs = [1, 2, 3]

colors = {
    "H3K27me3": "#027BB8",
    "H3K27ac": "#009697",
    "H3K4me3": "#E56F8E",
    "H3K9me3":"#a63c8d",
    "H3": "#808080"
}

for pc in pcs:
    # Extract PC1 and PC2
    pc1 = df_t[pc]
    pc2 = df_t[pc+1]

    for i, hue_color in zip([0,3],['day', 'histone']):

        if hue_color == 'histone':
            palette_color=colors
        else:
            palette_color='Set2'

        # Scatter plot: color by histone, shape by condition
        sns.scatterplot(x=pc1, y=pc2, hue=df_t.index.str.split('_').str[i] , style=df_t.index.str.split('_').str[2] , palette=palette_color)

        # Add labels and title
        plt.title(f'PC{pc} vs PC{pc+1} Scatter Plot')
        plt.xlabel(f'PC{pc}')
        plt.ylabel(f'PC{pc+1}')

        # Save the plot before showing it
        savepath = f"{file.removesuffix('.tsv')}_PC{pc}vsPC{pc+1}_{hue_color}_condition.pdf"
        print(f'fig saved to {savepath}')
        plt.savefig(savepath, bbox_inches='tight')

        # Show the plot
        plt.close()