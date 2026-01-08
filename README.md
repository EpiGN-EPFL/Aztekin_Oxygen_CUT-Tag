# Oxygen sensing & vertebrate limb regeneration

Tsissios et al., 2025

This repository contains analysis scripts for CUT&Tag experiments used in the study Oxygen sensing & vertebrate limb regeneration.

Reference genome mapping (i.e. generation of BAM and BigWig files) is performed using the DNAmapping pipelines from snakePipes. This repository focuses on downstream processing, QC, and visualization of the mapped data.

⸻

Repository structure
```
.
├── plotting
│   ├── bedmap.snakefile  # Quantify BigWig signals on predefined BED files
│   ├── bedmap_boxplot.py  # Plot bedmap results as boxplots
│   └── deeptools_plotProfile_plotHeatmap.snakefile  # deepTools plotProfile and plotHeatmap workflows
├── preprocessing  # Downstream of DNAmapping outputs
│   ├── bw_log2.snakefile  # Log2 transformation of BigWig files
│   ├── call_peaks_bam.snakefile  # Peak calling from BAM files
│   └── spikein_norm.snakefile  # Spike-in normalization
└── qc  # Additional QC (basic QC is included in DNAmapping)
    ├── frip.py  # Fragment-in-peaks (FRiP) score calculation
    ├── multiBWsummary_PCA.py  # Re-plot PCA from deepTools plotPCA outputs
    └── sampleQC_plot.py  # Helper functions for FRiP and spike-in QC plots

```
⸻

Dependencies

This repository assumes that alignment and initial processing have already been completed using snakePipes DNAmapping.

Core tools
- snakePipes (DNAmapping pipeline)
- Snakemake
- deepTools
- samtools
- bedtools

Python packages
- Python ≥ 3.8
- numpy
- pandas
- matplotlib
- seaborn

It is recommended to run the workflows in a conda environment consistent with your snakePipes setup.

⸻

Typical workflow
1.	Run DNAmapping (snakePipes)
Generate BAM and BigWig files using the DNAmapping pipeline.
2.	Preprocessing
- Perform spike-in normalization
- Log2-transform BigWig files
- Call peaks from BAM files
3.	Quality control
- Compute FRiP scores
- Inspect spike-in normalization
- Re-plot PCA from deepTools outputs
4.	Visualization
- Quantify signals on predefined BED regions
- Generate boxplots, profile plots, and heatmaps

⸻

Example usage

Run a Snakemake workflow:

```
snakemake -s preprocessing/bw_log2.snakefile --cores 8
```

Adjust paths and configuration files as needed for your dataset.

⸻

Citation

If you use this code, please cite:

Tsissios et al. (2025). 
