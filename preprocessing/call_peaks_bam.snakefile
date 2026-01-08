import glob
import os

samplenames, = glob_wildcards(os.path.join("filtered_bam",'{sample}.filtered.bam'))
peaks = ['broadPeak', 'narrowPeak']

rule all:
    input:
        expand("featureCounts/{sample}.{peak}_stats.txt", sample=samplenames, peak = peaks)

rule call_broad_peaks:
    input:
        bam="filtered_bam/{sample}.filtered.bam"
    output:
        broad_peaks=temp("peaks/{sample}_MACS3_bam.broadPeak.SAF")
    log:
        out="MACS3/.logs/call_broad_peaks.{sample}.out",
        err="MACS3/.logs/call_broad_peaks.{sample}.err"
    shell:
        """
        macs3 callpeak -t {input.bam} -f BAM --broad -g hs -n {wildcards.sample} \
        --outdir MACS3 -B --broad-cutoff 0.1 --nolambda >{log.out} 2>{log.err}
        awk 'OFS="\t" {{print $1"-"$2+1"-"$3, $1, $2+1, $3, "."}}' MACS3/{wildcards.sample}_peaks.broadPeak > {output.broad_peaks}
        """

rule call_narrow_peaks:
    input:
        bam="filtered_bam/{sample}.filtered.bam"
    output:
        narrow_peaks=temp("peaks/{sample}_MACS3_bam.narrowPeak.SAF")
    log:
        out="MACS3/.logs/call_narrow_peaks.{sample}.out",
        err="MACS3/.logs/call_narrow_peaks.{sample}.err"
    shell:
        """
        macs3 callpeak -t {input.bam} -f BAM -g hs -n {wildcards.sample} \
	    --outdir MACS3 -B -q 0.01 --nolambda >{log.out} 2>{log.err}
	    awk 'OFS="\t" {{print $1"-"$2+1"-"$3, $1, $2+1, $3, "."}}' MACS3/{wildcards.sample}_peaks.narrowPeak > {output.narrow_peaks}
        """

rule featureCounts:
    input:
        bam="filtered_bam/{sample}.filtered.bam",
        SAF="peaks/{sample}_MACS3_bam.{peak}.SAF"
    output:
        "featureCounts/{sample}.{peak}_stats.txt"
    log:
        "featureCounts/.logs/{sample}.{peak}.log"
    shell:
        """
        featureCounts -p --countReadPairs -a {input.SAF} -F SAF \
        -o {output} {input.bam} \
        --tmpDir /tmp/ \
        2> {log}
        """
        