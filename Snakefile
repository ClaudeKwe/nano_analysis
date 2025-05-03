# File: nano_analysis/Snakefile

# Configuration
configfile: "config_template.yaml"  # Default config, will be overridden by --configfile

import os
import glob

# Helper function to get samples from input directory structure
def get_samples(data_dir):
    samples = []
    for d in glob.glob(os.path.join(data_dir, "barcode*")):
        if os.path.isdir(d):
            sample = os.path.basename(d)
            samples.append(sample)
    return samples

# Get sample names from the input directory
SAMPLES = get_samples(config["input_dir"])

# Set results directory
RESULTS_DIR = config.get("results_dir", "results")

# Final output files that should be created
rule all:
    input:
        os.path.join(RESULTS_DIR, "mafft_alignment/all_samples_aligned.fasta"),
        os.path.join(RESULTS_DIR, "coverage_plots/all_samples_coverage.pdf"),
        expand(os.path.join(RESULTS_DIR, "consensus/{sample}.consensus.fasta"), sample=SAMPLES),
        os.path.join(RESULTS_DIR, "multiqc/multiqc_report.html")

# Combine all FASTQ files for each barcode/sample
rule combine_fastq:
    input:
        lambda wildcards: glob.glob(os.path.join(config["input_dir"], wildcards.sample, "*.fastq.gz"))
    output:
        temp(os.path.join(RESULTS_DIR, "combined_fastq/{sample}.fastq.gz"))
    shell:
        "cat {input} > {output}"

# Filter and trim reads using NanoFilt
rule nanofilt:
    input:
        os.path.join(RESULTS_DIR, "combined_fastq/{sample}.fastq.gz")
    output:
        os.path.join(RESULTS_DIR, "filtered_fastq/{sample}.filtered.fastq.gz")
    params:
        quality=config.get("min_read_quality", 8),
        length=config.get("min_read_length", 500),
        headcrop=config.get("headcrop", 50)
    conda:
        "envs/nanofilt.yaml"
    shell:
        "gunzip -c {input} | NanoFilt -q {params.quality} -l {params.length} --headcrop {params.headcrop} | gzip > {output}"

# Run FastQC on filtered reads
rule fastqc:
    input:
        fastq=os.path.join(RESULTS_DIR, "filtered_fastq/{sample}.filtered.fastq.gz")
    output:
        directory(os.path.join(RESULTS_DIR, "fastqc/{sample}"))
    conda:
        "envs/qc.yaml"
    shell:
        "mkdir -p {output} && fastqc {input.fastq} --outdir {output}"

# Align reads to reference using minimap2
rule minimap2_align:
    input:
        reads=os.path.join(RESULTS_DIR, "filtered_fastq/{sample}.filtered.fastq.gz"),
        ref=config["reference_genome"]
    output:
        os.path.join(RESULTS_DIR, "alignments/{sample}.bam")
    params:
        threads=config.get("threads", 8)
    conda:
        "envs/align.yaml"
    threads: lambda wildcards, params: params.threads
    shell:
        "minimap2 -ax map-ont -t {threads} {input.ref} {input.reads} | "
        "samtools sort -@ {threads} -o {output} - && "
        "samtools index {output}"

# Generate stats for the alignment
rule samtools_stats:
    input:
        os.path.join(RESULTS_DIR, "alignments/{sample}.bam")
    output:
        os.path.join(RESULTS_DIR, "stats/{sample}.stats.txt")
    conda:
        "envs/align.yaml"
    shell:
        "samtools stats {input} > {output}"

# Call variants using bcftools
rule call_variants:
    input:
        bam=os.path.join(RESULTS_DIR, "alignments/{sample}.bam"),
        ref=config["reference_genome"]
    output:
        vcf=os.path.join(RESULTS_DIR, "variants/{sample}.raw.vcf.gz")
    params:
        threads=config.get("threads", 4)
    conda:
        "envs/variants.yaml"
    threads: lambda wildcards, params: params.threads
    shell:
        "bcftools mpileup -Ou -f {input.ref} {input.bam} | "
        "bcftools call -mv -Oz -o {output.vcf} && "
        "bcftools index {output.vcf}"

# Filter variants for quality
rule filter_variants:
    input:
        os.path.join(RESULTS_DIR, "variants/{sample}.raw.vcf.gz")
    output:
        os.path.join(RESULTS_DIR, "variants/{sample}.filtered.vcf.gz")
    params:
        depth=config.get("variant_depth", 30),
        quality=config.get("variant_quality", 60)
    conda:
        "envs/variants.yaml"
    shell:
        "bcftools filter -i 'INFO/DP>{params.depth} && QUAL>{params.quality}' {input} -Oz -o {output} && "
        "bcftools index {output}"

# Generate bedgraph coverage file
rule generate_bedgraph:
    input:
        os.path.join(RESULTS_DIR, "alignments/{sample}.bam")
    output:
        os.path.join(RESULTS_DIR, "bedgraph/{sample}.bedgraph")
    conda:
        "envs/bedtools.yaml"
    shell:
        "genomeCoverageBed -ibam {input} -bg > {output}"

# Create coverage plots with R
rule plot_coverage:
    input:
        expand(os.path.join(RESULTS_DIR, "bedgraph/{sample}.bedgraph"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "coverage_plots/all_samples_coverage.pdf")
    params:
        results_dir=RESULTS_DIR
    conda:
        "envs/r_plots.yaml"
    script:
        "scripts/plot_coverage.R"

# Generate consensus sequence
rule consensus_sequence:
    input:
        bam=os.path.join(RESULTS_DIR, "alignments/{sample}.bam"),
        vcf=os.path.join(RESULTS_DIR, "variants/{sample}.filtered.vcf.gz"),
        ref=config["reference_genome"]
    output:
        temp=temp(os.path.join(RESULTS_DIR, "consensus/{sample}.temp.fasta")),
        final=os.path.join(RESULTS_DIR, "consensus/{sample}.consensus.fasta")
    conda:
        "envs/consensus.yaml"
    shell:
        """
        # Generate consensus
        bcftools consensus -f {input.ref} -o {output.temp} {input.vcf}
        
        # Replace header to include sample name
        awk '{{if(NR==1) print ">{wildcards.sample}"; else print $0}}' {output.temp} > {output.final}
        """

# Combine all consensus sequences
rule combine_consensus:
    input:
        ref=config["reference_genome"],
        consensus=expand(os.path.join(RESULTS_DIR, "consensus/{sample}.consensus.fasta"), sample=SAMPLES)
    output:
        temp(os.path.join(RESULTS_DIR, "mafft_alignment/all_consensus.fasta"))
    params:
        ref_name=config.get("reference_name", "reference")
    shell:
        """
        # Copy reference with proper name
        cat {input.ref} | awk '{{if($0 ~ /^>/) print ">{params.ref_name}"; else print $0}}' > {output}.ref
        
        # Combine with consensus sequences
        cat {output}.ref {input.consensus} > {output}
        rm {output}.ref
        """

# Align all consensus sequences with MAFFT
rule mafft_alignment:
    input:
        os.path.join(RESULTS_DIR, "mafft_alignment/all_consensus.fasta")
    output:
        os.path.join(RESULTS_DIR, "mafft_alignment/all_samples_aligned.fasta")
    params:
        threads=config.get("threads", 8)
    conda:
        "envs/mafft.yaml"
    threads: lambda wildcards,resources, input, attempt, params: params.threads
    shell:
        "mafft --thread {threads} --auto {input} > {output}"

# Run MultiQC to summarize all QC results
rule multiqc:
    input:
        expand(os.path.join(RESULTS_DIR, "fastqc/{sample}"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "stats/{sample}.stats.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "multiqc/multiqc_report.html")
    params:
        fastqc_dir=os.path.join(RESULTS_DIR, "fastqc"),
        stats_dir=os.path.join(RESULTS_DIR, "stats"),
        out_dir=os.path.join(RESULTS_DIR, "multiqc")
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc {params.fastqc_dir} {params.stats_dir} -o {params.out_dir}"
