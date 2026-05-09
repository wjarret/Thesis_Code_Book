#RNP_MaP Snakemake Pipelne
#Created February 2026
#@author: wjarret and chat GPT 5.2

configfile: "config_RNP_MaP.yaml"

fastp = config["fastp"]
samples_config = config.get("samples", [])
Index_adapters = config["Index_adapters"]
Cutadapt = config["Cutadapt"]
seqkit = config["seqkit"]


# Discover samples from filesystem
 

(samples,) = glob_wildcards("data/samples/{sample}_R1_001.fastq.gz")

#print(f"Total samples discovered: {len(samples)}")
#print(f"All samples: {samples}")

DMSO_SAMPLES = [s for s in samples if "DMSOaq" in s]
OOPS_SAMPLES = [s for s in samples if "OOPS" in s]

#print(f"DMSO samples for variant calling: {len(DMSO_SAMPLES)}")
#print(f"OOPS samples for preprocessing only: {len(OOPS_SAMPLES)}")
#print(f"DMSO samples: {DMSO_SAMPLES}")
#print(f"OOPS samples: {OOPS_SAMPLES}")


 
# Final target
 

rule all:
    input:
        # DMSO runs full pipeline
        expand("primer_reinsertion/primer_reinsertion_ASV_{sample}.fa", sample=DMSO_SAMPLES),
        # OOPS stops at UMI extraction
        expand("umi_removed/{sample}_R1_001.fastq.gz", sample=OOPS_SAMPLES),
        expand("umi_removed/{sample}_R2_001.fastq.gz", sample=OOPS_SAMPLES)

 
# Copy FASTQ
 

rule copy_fastq:
    input:
        r1="data/samples/{sample}_R1_001.fastq.gz",
        r2="data/samples/{sample}_R2_001.fastq.gz"
    output:
        r1="raw/{sample}_R1_001.fastq.gz",
        r2="raw/{sample}_R2_001.fastq.gz"
    shell:
        """
        mkdir -p raw logs
        cp {input.r1} {output.r1}
        cp {input.r2} {output.r2}
        """

 
# Adapter trimming step 1
 

rule cutadapt_step1:
    input:
        r1="raw/{sample}_R1_001.fastq.gz",
        r2="raw/{sample}_R2_001.fastq.gz"
    output:
        r1="index_adapter_trimming/step1_{sample}_R1_001.fastq.gz",
        r2="index_adapter_trimming/step1_{sample}_R2_001.fastq.gz"
    threads: 8
    params:
        adapter_NT=Index_adapters["Non_Transposase"],
        error=Cutadapt["error_rate"]
    log:
        "logs/cutadapt_step1/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p index_adapter_trimming logs/cutadapt_step1
        cutadapt \
          -j {threads} \
          -a file:{params.adapter_NT} \
          -A file:{params.adapter_NT} \
          -e {params.error} \
          --report=minimal \
          -o {output.r1} -p {output.r2} \
          {input.r1} {input.r2} > {log} 2>&1
        """

 
# Adapter trimming step 2
 

rule cutadapt_step2:
    input:
        r1="index_adapter_trimming/step1_{sample}_R1_001.fastq.gz",
        r2="index_adapter_trimming/step1_{sample}_R2_001.fastq.gz"
    output:
        r1="index_adapter_trimming/step2_{sample}_R1_001.fastq.gz",
        r2="index_adapter_trimming/step2_{sample}_R2_001.fastq.gz"
    threads: 8
    params:
        adapter_T=Index_adapters["Transposase"],
        error=Cutadapt["error_rate"]
    log:
        "logs/cutadapt_step2/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p index_adapter_trimming logs/cutadapt_step2
        cutadapt \
            -j {threads} \
            -a file:{params.adapter_T} \
            -A file:{params.adapter_T} \
            -e {params.error} \
            --report=minimal \
            -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} > {log} 2>&1
        """

 
# UMI extraction
 

rule umi_extraction:
    input:
        r1="index_adapter_trimming/step2_{sample}_R1_001.fastq.gz",
        r2="index_adapter_trimming/step2_{sample}_R2_001.fastq.gz"
    output:
        r1="umi_removed/{sample}_R1_001.fastq.gz",
        r2="umi_removed/{sample}_R2_001.fastq.gz"
    threads: 2
    params:
        umi_len=fastp["umi_len"]
    log:
        "logs/fastp_umi/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p umi_removed logs/fastp_umi
        fastp \
            --thread {threads} \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --umi \
            --umi_loc per_read \
            --umi_len {params.umi_len} \
            --umi_prefix UMI \
            -j umi_removed/{wildcards.sample}.json \
            -h umi_removed/{wildcards.sample}.html \
            > {log} 2>&1
        """

 
# Organize DMSO samples
 

rule organize_DMSO:
    input:
        r1="umi_removed/{sample}_R1_001.fastq.gz",
        r2="umi_removed/{sample}_R2_001.fastq.gz"
    output:
        r1="conditioned/DMSO/{sample}_R1_001.fastq.gz",
        r2="conditioned/DMSO/{sample}_R2_001.fastq.gz"
    shell:
        """
        mkdir -p conditioned/DMSO
        ln -sf $(realpath {input.r1}) {output.r1}
        ln -sf $(realpath {input.r2}) {output.r2}
        """

 
# Primer trimming
 

rule cutadapt_primer_removal:
    input:
        r1="conditioned/DMSO/{sample}_R1_001.fastq.gz",
        r2="conditioned/DMSO/{sample}_R2_001.fastq.gz"
    output:
        r1="cutadapt_primer_removal/{sample}_R1_001.fastq.gz",
        r2="cutadapt_primer_removal/{sample}_R2_001.fastq.gz"
    threads: 4
    params:
        primer_1=Cutadapt["R1_primer_Length"],
        primer_2=Cutadapt["R2_primer_Length"],
        min_length=Cutadapt["min_length"],
        max_length=Cutadapt["max_length"]
    log:
        "logs/cutadapt_primer_removal/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p cutadapt_primer_removal logs/cutadapt_primer_removal
        cutadapt \
            -j {threads} \
            -u {params.primer_1} \
            -U {params.primer_2} \
            --minimum-length {params.min_length} \
            --maximum-length {params.max_length} \
            --report=minimal \
            -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} > {log} 2>&1
        """

 
# BBMerge
 

rule BB_Merge:
    input:
        r1="cutadapt_primer_removal/{sample}_R1_001.fastq.gz",
        r2="cutadapt_primer_removal/{sample}_R2_001.fastq.gz"
    output:
        r1="BB_Merge/{sample}_R1_001.fastq.gz",
        r2="BB_Merge/{sample}_R2_001.fastq.gz"
    log:
        "logs/BB_Merge/{sample}.log"
    shell:
        """
        mkdir -p BB_Merge logs/BB_Merge
        /home/wjarret/micromamba/bbmap/bbmerge.sh \
            in1={input.r1} \
            in2={input.r2} \
            out1={output.r1} \
            out2={output.r2} \
            tbo=t \
            merge=false > {log} 2>&1
        """

 
# Quality trimming
 

rule Quality_trimming_fastp:
    input:
        r1="BB_Merge/{sample}_R1_001.fastq.gz",
        r2="BB_Merge/{sample}_R2_001.fastq.gz"
    output:
        r1="Quality_trimming_fastp/{sample}_R1_001.fastq.gz",
        r2="Quality_trimming_fastp/{sample}_R2_001.fastq.gz"
    threads: 4
    params:
        quality_score=fastp["quality_score_threshold"]
    log:
        "logs/Quality_trimming_fastp/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p Quality_trimming_fastp logs/Quality_trimming_fastp
        fastp \
          -i {input.r1} \
          -I {input.r2} \
          -o {output.r1} \
          -O {output.r2} \
          --detect_adapter_for_pe \
          --correction \
          -q {params.quality_score} \
          --n_base_limit 0 \
          --thread {threads} > {log} 2>&1
        """

 
# DADA2
 

rule dada2:
    input:
        r1="Quality_trimming_fastp/{sample}_R1_001.fastq.gz",
        r2="Quality_trimming_fastp/{sample}_R2_001.fastq.gz"
    output:
        fa="dada2/ASV_{sample}.fa"
    log:
        "logs/dada2/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p dada2 logs/dada2
        Rscript scripts/dada2_snakemake.R \
            --r1 {input.r1} \
            --r2 {input.r2} \
            --out {output.fa} > {log} 2>&1
        """

 
# Primer reinsertion
 

rule primer_reinsertion:
    input:
       "dada2/ASV_{sample}.fa"
    output:
        "primer_reinsertion/primer_reinsertion_ASV_{sample}.fa"
    params:
        adapter_f=seqkit["adapter_forward"]
    log:
        "logs/primer_reinsertion/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        """
        mkdir -p primer_reinsertion logs/primer_reinsertion
        seqkit mutate --insertion 0:{params.adapter_f} {input} -o {output} > {log} 2>&1
        """
