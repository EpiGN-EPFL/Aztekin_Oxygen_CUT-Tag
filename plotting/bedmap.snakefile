import os
import pandas as pd
from snakemake.io import glob_wildcards

# -----------------------------------------------------------------------------
# CONFIG & GLOBAL VARIABLES
# -----------------------------------------------------------------------------

# Directory containing the bigWig files (set in your config.yaml)
indir = config["indir"]

# Discover all region BED files (without extension)
# e.g. regions_bedmap/myregion.bed -> bedfiles = ["myregion", ...]
bedfiles, = glob_wildcards(os.path.join("regions_bedmap", "{bed}.bed"))
print("Regions:", bedfiles)

# Discover all sample bigWig files in `indir`
# e.g. indir/sample1.bw -> samples = ["sample1", ...]
samples, = glob_wildcards(os.path.join(indir, "{sample}.bw"))
print("Samples:", samples)

# Derived directories (helps keep things consistent)
BED_FROM_BW_DIR      = f"bed_{indir}"
BED_FROM_BW_LOG_DIR  = os.path.join(BED_FROM_BW_DIR, "logs")
BEDMAP_DIR           = f"bedmap_{indir}"
BEDMAP_LOG_DIR       = os.path.join(BEDMAP_DIR, "logs")
COMBINED_BEDMAP_DIR  = f"combined_bedmap_{indir}"

# -----------------------------------------------------------------------------
# RULE: all
# -----------------------------------------------------------------------------
rule all:
    input:
        # All sample x region bedmap outputs
        expand(
            os.path.join(BEDMAP_DIR, "{sample}_on_{region}.bed"),
            sample=samples,
            region=bedfiles,
        ),
        # Combined bedmaps per region
        expand(
            os.path.join(COMBINED_BEDMAP_DIR, "{region}_combined.bed"),
            region=bedfiles,
        )

# -----------------------------------------------------------------------------
# RULE: bw_to_bed
#   Convert bigWig signal to BED using bigWigToBedGraph.
# -----------------------------------------------------------------------------
rule bw_to_bed:
    input:
        B_bw = os.path.join(indir, "{sample}.bw")
    output:
        B_bed = os.path.join(BED_FROM_BW_DIR, "{sample}.bed")
    log:
        os.path.join(BED_FROM_BW_LOG_DIR, "{sample}_BWtoBed.log")
    shell:
        r"""
        mkdir -p {os.path.dirname(output.B_bed)} {os.path.dirname(log)}
        bigWigToBedGraph {input.B_bw} {output.B_bed} 2> {log}
        """

# -----------------------------------------------------------------------------
# RULE: bed_map
#   Map per-base/interval signal from sample BED onto region BED via bedtools map.
#
#   The output has columns:
#     chr, start, end, {sample}_sum, {sample}_mean, {sample}_median
# -----------------------------------------------------------------------------
rule bed_map:
    input:
        A_bed = "regions_bedmap/{region}.bed",
        B_bed = os.path.join(BED_FROM_BW_DIR, "{sample}.bed")
    output:
        bedmap = os.path.join(BEDMAP_DIR, "{sample}_on_{region}.bed")
    log:
        os.path.join(BEDMAP_LOG_DIR, "{sample}_{region}.log")
    params:
        # bedtools map options: summarize column 4 of B_bed with sum/mean/median
        options = "-c 4 -o sum,mean,median"
    shell:
        r"""
        mkdir -p {os.path.dirname(output.bedmap)} {os.path.dirname(log)}

        # Write header with sample-specific column names
        echo -e "chr\tstart\tend\t{wildcards.sample}_sum\t{wildcards.sample}_mean\t{wildcards.sample}_median" > {output.bedmap}

        # Pipe regions (chr, start, end) into bedtools map
        awk '{{print $1, $2, $3}}' OFS="\t" {input.A_bed} \
        | bedtools map -a - -b {input.B_bed} {params.options} >> {output.bedmap} 2> {log}
        """

# -----------------------------------------------------------------------------
# RULE: combine_bedmaps
#   Combine all sample-specific bedmap files for a given region into one table.
# -----------------------------------------------------------------------------
rule combine_bedmaps:
    input:
        # All bedmaps corresponding to the same region across samples
        lambda wildcards: expand(
            os.path.join(BEDMAP_DIR, "{sample}_on_{region}.bed"),
            sample=samples,
            region=wildcards.region,
        )
    output:
        os.path.join(COMBINED_BEDMAP_DIR, "{region}_combined.bed")
    run:
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        # Read each bedmap file into a DataFrame
        bedmap_dfs = [pd.read_csv(bf, sep="\t") for bf in input]

        # Iteratively merge on the genomic coordinates
        merged_df = bedmap_dfs[0]
        for df in bedmap_dfs[1:]:
            merged_df = pd.merge(
                merged_df,
                df,
                on=["chr", "start", "end"],
                how="outer",
            )

        # Write combined table
        merged_df.to_csv(output[0], sep="\t", index=False)