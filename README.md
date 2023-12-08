# IRE1aKO_Polysome_RNAseq

This is the repository for analysis of bulk RNAseq from P0 mouse cerebral cortex, isolated either as 
- 1. Total RNAseq (total RNA library with rRNA depletion).
- 2. Polysomes (enrichment for polysomes, followed by RNA extraction and total RNA library with rRNA depletion)
- 3. Heavy Polysomes (enrichment for polysomes, selected 'heavy' fraction, followed by RNA extraction and total RNA library with rRNA depletion).
- 4. light Polysomes (enrichment for polysomes, selected 'light' fraction, followed by RNA extraction and total RNA library with rRNA depletion).

These isolation and their batches are provided in the 'metadata' object described in the main `Analysis.Rmd` script, which contains all code written for statistical analysis and plotting.


# Mapping and Quantification

Fastq files specified in `IRE1a_Poly_config.yaml` were quantified by Salmon and additionally aligned to GRCm39 using STAR using the following snakefile:

```{python Snakefile}
#! bin/env/python


configfile: "/fast/users/newmana_c/config/IRE1a_Poly_config.yaml"

TE_GTFFILE = config["TE_gtf_file"]
GTFFILE = config["gene_gtf_file"]
STAR_INDEX = config["star_index"]


import pandas as pd

# Flatten nested dictionary and create a dataframe
sample_data = pd.DataFrame([
    {"sample": sample, "read_type": rt, **details}
    for rt, samples in config["samples"].items()
    for sample, details in samples.items()
])

SAMPLES_PE = sample_data.query("read_type == 'PE'")["sample"].tolist()
SAMPLES_SE = sample_data.query("read_type == 'SE'")["sample"].tolist()

# This function will look up `name` for a given `sample` in the dataframe
def get_name(wildcards):
    return sample_data.query(f"sample == '{wildcards.sample}'")["name"].values[0]

SE_INDEX =expand("alignment/SE/{sample}.sorted.bam.bai", sample = SAMPLES_SE),
PE_INDEX =expand("alignment/PE/{sample}.sorted.bam.bai", sample = SAMPLES_PE),

rule all:
    input:
        SE_INDEX =expand("alignment/SE/{sample}.sorted.bam.bai", sample = SAMPLES_SE),
        PE_INDEX =expand("alignment/PE/{sample}.sorted.bam.bai", sample = SAMPLES_PE),
        PE_SWIM = expand("salmon_quant/PE/{sample}", sample=SAMPLES_PE),
        SE_SWIM = expand("salmon_quant/SE/{sample}", sample=SAMPLES_SE),
        SE_TET_COUNTS = expand("counts/SE/{sample}.cntTable", sample = SAMPLES_SE),
        PE_TET_COUNTS = expand("counts/PE/{sample}.cntTable", sample = SAMPLES_PE)




rule swim_PE:
    input:
        fq1 = lambda wildcards: f"fastq/concatenated/P1757_RNA_{get_name(wildcards)}_L001_R1_001.fastq.gz",
        fq2 = lambda wildcards: f"fastq/concatenated/P1757_RNA_{get_name(wildcards)}_L001_R2_001.fastq.gz",
        index = "/fast/users/newmana_c/work/genomes/GRCm39/salmon1.4_gencode.vM33"
    output:
        directory("salmon_quant/PE/{sample}")
    log:
        "logs/PE/{sample}_align_sort.log"
    threads: 8
    resources:
        mem_mb=24000
    params:
        lib_type = 'A'
    shell:
        """
        salmon quant -i {input.index} \
              -p {threads} \
              --validateMappings \
              -l {params.lib_type} \
              -1 {input.fq1} \
              -2 {input.fq2} \
              -o {output}
        """

rule swim_SE:
    input:
        fq = lambda wildcards: f"fastq/concatenated/P1757_RNA_{get_name(wildcards)}_L001_R1_001.fastq.gz",
        index = "/fast/users/newmana_c/work/genomes/GRCm39/salmon1.4_gencode.vM33"
    output:
        directory("salmon_quant/SE/{sample}")
    log:
        "logs/SE/{sample}_align_sort.log"
    threads: 8
    resources:
        mem_mb=24000
    params:
        lib_type = 'A'
    shell:
        """
        salmon quant -i {input.index} \
              -p {threads} \
              -l {params.lib_type} \
              --validateMappings \
              -r {input.fq} \
              -o {output}
        """


rule index_bam_PE:
	input:
		"alignment/PE/{sample}.sorted.bam"
	output:
		"alignment/PE/{sample}.sorted.bam.bai"
	log:
                "logs/{sample}.index_bam"
	threads: 1
        resources:
            mem='10G',
            time='04:00:00'
	shell:
                """
                samtools index {input} 2> {log}
                """

rule index_bam_SE:
	input:
		"alignment/SE/{sample}.sorted.bam"
	output:
		"alignment/SE/{sample}.sorted.bam.bai"
	log:
                "logs/{sample}.index_bam"
	threads: 1
        resources:
            mem='10G',
            time='04:00:00'
	shell:
                """
                samtools index {input} 2> {log}
                """

rule TET_pe_counts:
        input:
                gtf = GTFFILE,
                te = TE_GTFFILE,
                BAI = PE_INDEX,
                bam = "alignment/PE/{sample}.sorted.bam"
        output: "counts/PE/{sample}.cntTable"
        threads: 8
        resources:
            mem_mb=24000
        shell:"""
              TEcount \
              --mode multi \
              --format BAM \
              --sortByPos \
              -b {input.bam} \
              --GTF {input.gtf} \
              --TE {input.te} \
              --project counts/PE/{wildcards.sample}
              """





rule TET_se_counts:
        input:
                gtf = GTFFILE,
                te = TE_GTFFILE,
                BAI = SE_INDEX,
                bam = "alignment/SE/{sample}.sorted.bam"
        output: "counts/SE/{sample}.cntTable"
        threads: 8
        resources:
            mem_mb=24000
        shell:"""
              TEcount \
              --mode multi \
              --format BAM \
              --sortByPos \
              -b {input.bam} \
              --GTF {input.gtf} \
              --TE {input.te} \
              --project counts/SE/{wildcards.sample}
              """

rule sort_PE:
        input: "alignment/PE/{sample}_Aligned.out.bam"
        output: "alignment/PE/{sample}.sorted.bam"
        threads: 4
        resources: 
            mem='5G',
            time='30:00'
        shell:"""
                samtools sort -m 1G -@ {threads} -O bam -T {output}.tmp {input} -o {output};
                rm alignment/PE/{wildcards.sample}_Aligned.out.bam
                """

rule sort_SE:
        input: "alignment/SE/{sample}_Aligned.out.bam"
        output: "alignment/SE/{sample}.sorted.bam"
        threads: 4
        resources: 
            mem='5G',
            time='30:00'
        shell:"""
                samtools sort -m 1G -@ {threads} -O bam -T {output}.tmp {input} -o {output};
                rm alignment/SE/{wildcards.sample}_Aligned.out.bam
                """



rule align_PE:
   input:
       fq1=lambda wildcards: f"fastq/concatenated/P1757_RNA_{get_name(wildcards)}_L001_R1_001.fastq.gz",
       fq2=lambda wildcards: f"fastq/concatenated/P1757_RNA_{get_name(wildcards)}_L001_R2_001.fastq.gz",
       genome=STAR_INDEX,
       gtf=GTFFILE
   output:"alignment/PE/{sample}_Aligned.out.bam"
   log: "logs/{sample}_align_sort.log"
   threads: 16
   resources:
       mem_mb=36384,
       time='04:00:00'
   shell:"""
       STAR --genomeDir {input.genome} \
       --outFileNamePrefix alignment/PE/{wildcards.sample}_ \
       --readFilesIn {input.fq1} {input.fq2} \
       --readFilesCommand zcat \
       --runThreadN {threads} \
       --genomeLoad NoSharedMemory \
       --outSAMattributes All \
       --outFilterMultimapNmax 100 \
       --winAnchorMultimapNmax 100 \
       --outSAMstrandField intronMotif \
       --outSAMtype BAM Unsorted \
       --sjdbGTFfile {input.gtf} 
       mkdir -p starlogs
       mv -f alignment/PE/{wildcards.sample}_Log.final.out alignment/PE/{wildcards.sample}_Log.out alignment/PE/{wildcards.sample}_Log.progress.out alignment/PE/{wildcards.sample}_SJ.out.tab alignment/PE/{wildcards.sample}__STARgenome starlogs
       """


rule align_SE:
   input:
       fq1=lambda wildcards: f"fastq/concatenated/P1757_RNA_{get_name(wildcards)}_L001_R1_001.fastq.gz",
       genome=STAR_INDEX,
       gtf=GTFFILE
   output:"alignment/SE/{sample}_Aligned.out.bam"
   log: "logs/{sample}_align_sort.log"
   threads: 16
   resources:
       mem_mb=36384,
       time='04:00:00'
   shell:"""
       STAR --genomeDir {input.genome} \
       --outFileNamePrefix alignment/SE/{wildcards.sample}_ \
       --readFilesIn {input.fq1} \
       --readFilesCommand zcat \
       --runThreadN {threads} \
       --genomeLoad NoSharedMemory \
       --outSAMattributes All \
       --outFilterMultimapNmax 100 \
       --winAnchorMultimapNmax 100 \
       --outSAMstrandField intronMotif \
       --outSAMtype BAM Unsorted \
       --sjdbGTFfile {input.gtf} 
       mkdir -p starlogs
       mv -f alignment/SE/{wildcards.sample}_Log.final.out alignment/SE/{wildcards.sample}_Log.out alignment/SE/{wildcards.sample}_Log.progress.out alignment/SE/{wildcards.sample}_SJ.out.tab alignment/SE/{wildcards.sample}__STARgenome starlogs
       """

```

