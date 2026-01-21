# Updated Drosophila Variant Discovery Pipeline
SAMPLES = ["sample"]

rule all:
    input:
        # Target the FILTERED VCF now
        expand("results/vcf/{sample}_filtered.vcf", sample=SAMPLES),
        "results/qc/multiqc_report.html"

rule align:
    input:
        ref="data/ref/genome.fa",
        r1="data/raw/{sample}_R1.fastq",
        r2="data/raw/{sample}_R2.fastq"
    output:
        bam="results/vcf/{sample}.bam"
    threads: 16
    shell:
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | "
        "samtools view -Sb - | samtools sort -o {output.bam}"

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
        vcf="results/vcf/{sample}_raw.vcf" # Renamed to 'raw'
    shell:
        "freebayes -f {input.ref} {input.bam} > {output.vcf}"

# --- NEW: Phase 9 Filtering Rule ---
rule filter_variants:
    input:
        vcf="results/vcf/{sample}_raw.vcf"
    output:
        filtered="results/vcf/{sample}_filtered.vcf"
    shell:
        # Standard filters for Drosophila: 
        # QUAL > 30 (99.9% accuracy) and DP > 10 (sufficient coverage)
        "vcffilter -f 'QUAL > 30 & DP > 10' {input.vcf} > {output.filtered}"

rule multiqc:
    input:
        expand("results/vcf/{sample}_filtered.vcf", sample=SAMPLES)
    output:
        "results/qc/multiqc_report.html"
    shell:
        "multiqc results/ -o results/qc/ --force"