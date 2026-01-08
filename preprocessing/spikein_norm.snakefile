import glob, os

configfile: os.path.join("DNAmapping_organism.yaml") # config
spikein_name='NC_001416.1' #### maybe change this to config
samplenames, = glob_wildcards(os.path.join("filtered_bam",'{sample}.filtered.bam'))

rule all:
    input:
        os.path.join("bamCoverage_scaleFactor","multiBamSummary.spike_in.scaleFactors.tsv"),
        expand(os.path.join("bamCoverage_scaleFactor","{sample}.mapq{mapq}.scaled.bw"), sample = samplenames, mapq=[3])


rule spike_in_region:
    input:
        fai=config['genome_index']
    params:
        spikein=spikein_name
    output:
        bed=os.path.join("bamCoverage_scaleFactor","spike_in.bed")
    shell:
        "grep '{params.spikein}' {input.fai} | awk -v OFS='\t' '{{print($1,0,$2)}}' > {output.bed}"

rule scalingFactors:
    input:
        bam_files=expand(os.path.join("filtered_bam","{sample}.filtered.bam"), sample = samplenames),
        bai_files=expand(os.path.join("filtered_bam","{sample}.filtered.bam.bai"), sample = samplenames),
        spikein_bed=rules.spike_in_region.output.bed
    output:
        outtable=os.path.join("bamCoverage_scaleFactor","multiBamSummary.spike_in.npz"),
        rawtable=os.path.join("bamCoverage_scaleFactor","multiBamSummary.spike_in.rawcounts.tsv"),
        factortable=os.path.join("bamCoverage_scaleFactor","multiBamSummary.spike_in.scaleFactors.tsv")
    params:
        binsize="--binSize 1000"
    threads: 32
    log: os.path.join("bamCoverage_scaleFactor", "log", "multiBamSummary.bins.log")
    shell:
        "multiBamSummary bins -p {threads} {params.binsize} --region $(cat {input.spikein_bed} | tr '\t' ':') --bamfiles {input.bam_files} --outFileName {output.outtable} --scalingFactors {output.factortable} --outRawCounts {output.rawtable} &> {log}"

# rule bamCoverage_scaleFactor:
#     input:
#         bam=os.path.join("filtered_bam","{sample}.filtered.bam"),
#         factortable=rules.scalingFactors.output.factortable
#     output:
#         bigwig=os.path.join("bamCoverage_scaleFactor", "{sample}.scaled.bw")
#     params:
#         param_string="--binSize 25 --extendReads",
#         pattern="{sample}.filtered.bam"
#     threads: 8
#     log: os.path.join("bamCoverage_scaleFactor", "log", "bamCoverage.{sample}.scaled.log")
#     shell:
#         "bamCoverage -p {threads} {params.param_string} -b {input.bam} -o {output.bigwig} --scaleFactor $(grep '{params.pattern}' {input.factortable} | cut -f 2) &> {log}"

rule bamCoverage_scaleFactor:
    input:
        bam=os.path.join("filtered_bam","{sample}.filtered.bam"),
        factortable=rules.scalingFactors.output.factortable
    output:
        bigwig=os.path.join("bamCoverage_scaleFactor", "{sample}.mapq{mapq}.scaled.bw")
    params:
        param_string="--binSize 25 --minMappingQuality {mapq} --extendReads",
        pattern="{sample}.filtered.bam"
    threads: 8
    log: os.path.join("bamCoverage_scaleFactor", "log", "bamCoverage.{sample}.mapq{mapq}.scaled.log")
    shell:
        "bamCoverage -p {threads} {params.param_string} -b {input.bam} -o {output.bigwig} --scaleFactor $(grep '{params.pattern}' {input.factortable} | cut -f 2) &> {log}"
