import glob
import os

# Read config
indir=config['indir']
genome_size=config['genome']
outdir=f'{indir}_log2'

# get all samples
samplenames, = glob_wildcards(os.path.join(indir,'{sample}.bw'))
print(samplenames)


rule all:
    input:
        expand("{dir}_log2/{sample}.log2.bw", dir=indir, sample=samplenames)

rule log2_bed:
    input:
        indir+"/{sample}.bw"
    output:
        "bed_"+outdir+"/{sample}.log2.bed",
    log:
        outdir+"/logs/{sample}.log2bed.log"
    shell:
        """
        wiggletools offset 1 {input} | wiggletools write_bg {output} log 2 - 2> {log}
        """

rule bed_to_bw:
    input:
        "bed_"+outdir+"/{sample}.log2.bed"
    output:
        outdir+"/{sample}.log2.bw"
    params:
        genome_size=genome_size,
    log:
        outdir+"/logs/{sample}.log2bw.log"
    shell:
        """
        bedGraphToBigWig {input} {params.genome_size} {output} 2>{log}
        """