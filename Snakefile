SAMPLES = ["sample"]

rule all:
    input:
        expand("results/vcf/{sample}_filtered.vcf", sample=SAMPLES),
        "results/qc/multiqc_report.html"

rule align:
    input:
        ref="data/ref/genome.fa",
        r1="data/raw/{sample}_R1.fastq",
        r2="data/raw/{sample}_R2.fastq"
    output:
        bam="results/vcf/{sample}.bam"
    log:
        "results/qc/{sample}_bwa_mem.log"
    threads: 16
    shell:
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} 2> {log} | "
        "samtools view -Sb - | samtools sort -o {output.bam}"

rule bam_stats:
    input:
        bam="results/vcf/{sample}.bam"
    output:
        flagstat="results/qc/{sample}.flagstat",
        stats="results/qc/{sample}.stats"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat} && "
        "samtools stats {input.bam} > {output.stats}"

rule index_bam:
    input: "results/vcf/{sample}.bam"
    output: "results/vcf/{sample}.bam.bai"
    shell: "samtools index {input}"

rule call_variants:
    input:
        ref="data/ref/genome.fa",
        bam="results/vcf/{sample}.bam",
        bai="results/vcf/{sample}.bam.bai"
    output:
        vcf="results/vcf/{sample}_raw.vcf"
    log:
        "results/qc/{sample}_freebayes.log"
    shell:
        "freebayes -f {input.ref} {input.bam} > {output.vcf} 2> {log}"

rule filter_variants:
    input:
        vcf="results/vcf/{sample}_raw.vcf"
    output:
        filtered="results/vcf/{sample}_filtered.vcf"
    shell:
        "vcffilter -f 'QUAL > 30 & DP > 10' {input.vcf} > {output.filtered}"

rule multiqc:
    input:
        vcf=expand("results/vcf/{sample}_filtered.vcf", sample=SAMPLES),
        bam_stats=expand("results/qc/{sample}.stats", sample=SAMPLES),
        flagstat=expand("results/qc/{sample}.flagstat", sample=SAMPLES),
        bwa_log=expand("results/qc/{sample}_bwa_mem.log", sample=SAMPLES),
        fb_log=expand("results/qc/{sample}_freebayes.log", sample=SAMPLES)
    output:
        "results/qc/multiqc_report.html"
    shell:
        "multiqc results/qc/ -o results/qc/ --force"