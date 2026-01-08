import glob
import os

colors = {
    "H3K27me3": "#027BB8",
    "H3K27ac": "#009697",
    "H3K4me3": "#E56F8E",
    "H3K9me3":"#a63c8d",
    "H3": "#050505"
}
indir=config['indir']
outdir=f'profile_{indir}'
histones=config['histones']
e=config['e']
print(histones)
bed_files, = glob_wildcards(os.path.join("regions_profile", "{bed}.bed"))
print(bed_files)

sort_priority = {
    "ALI": 0,
    "96-well": 1
}

def custom_sort_key(filename):
    # Determine the primary sort order
    for key in sort_priority:
        if key in filename:
            return (sort_priority[key], filename)
    # Return a large value and the filename as fallback
    return (len(sort_priority), filename)


# for plotHeatmap, sort the bw first
def get_bw(indir, histones, e):
    bws = list()
    color_list = list()
    for hist in histones:
        bw = [f for f in glob.glob(f'{indir}/*{hist}{e}*.bw') if 'H3_sub.H3' not in f]
        color_list += [f"'#ffffff,{colors[hist]}'"] * len(bw)  # Create alternating color pairs
        bw = sorted(bw, key=custom_sort_key)
        bws += bw
    color_string = " ".join(color_list)  # Join the color list into a single string
    print(color_string)
    return bws, color_string  # Return the list of files and the color string


rule all:
    input:
        expand("{dir}/{bed}_{cmd}_profile.pdf", dir=outdir, bed=bed_files, cmd=['referencePoint','scaleRegions']),
        expand("{dir}/{bed}_{cmd}_heatmap.pdf", dir=outdir, bed=bed_files, cmd=['referencePoint','scaleRegions'])


rule referencePoint_computeMartix:
    input:
        bedGraph="regions_profile/{bed}.bed",
        bigWig=get_bw(indir, histones, e)[0]
    output:
        "{dir}/{bed}_referencePoint_mat.gz",
    log:
        "{dir}/logs/{bed}_referencePoint_computeMatrix.log"
    threads: 8
    params:
        regions="-b 2000 -a 2000",
        bins="--binSize 50"
    shell:
        """
        computeMatrix reference-point \
        --referencePoint center \
        -R {input.bedGraph} \
        -S {input.bigWig} \
        --missingDataAsZero \
        --skipZeros \
        -o {output} \
        --smartLabels \
        -p max/2 \
        {params.bins} {params.regions} 2>{log}
        """

rule scaleRegions_computeMartix:
    input:
        bedGraph="regions_profile/{bed}.bed",
        bigWig=get_bw(indir, histones, e)[0]
    output:
        "{dir}/{bed}_scaleRegions_mat.gz",
    log:
        "{dir}/logs/{bed}_scaleRegions_computeMatrix.log"
    threads: 8
    params:
        regions="-b 5000 -a 5000 --regionBodyLength 2000",
        bins="--binSize 100"
    shell:
        """
        computeMatrix scale-regions \
        -R {input.bedGraph} \
        -S {input.bigWig} \
        --missingDataAsZero \
        --skipZeros \
        -o {output} \
        --smartLabels \
        -p max/2 \
        {params.bins} {params.regions} 2>{log}
        """


rule plotProfile:
    input:
        matrix="{dir}/{bed}_{cmd}_mat.gz"
    output:
        fig="{dir}/{bed}_{cmd}_profile.pdf",
        bed="{dir}/{bed}_{cmd}_profile.bed",
    log:
        "{dir}/logs/{bed}_{cmd}_plotProfile.log"
    shell:
        """
        plotProfile -m {input.matrix} \
        -o {output.fig} \
        --refPointLabel "peak" \
        --outFileSortedRegions {output.bed} 2>{log}
        """

rule plotHeatmap:
    input:
        matrix="{dir}/{bed}_{cmd}_mat.gz"
    output:
        fig="{dir}/{bed}_{cmd}_heatmap.pdf",
    log:
        "{dir}/logs/{bed}_{cmd}_plotHeatmap.log"
    params:
        colorscale="--zMin 0 --zMax 4",
        color=get_bw(indir, histones, e)[1]
    shell:
        """
        plotHeatmap -m {input.matrix} \
        -o {output.fig} \
        --refPointLabel "peak" \
        --colorList {params.color} 2>{log}
        """