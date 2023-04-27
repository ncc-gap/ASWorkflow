configfile: "config.yaml"

chr_nums = config["chromosomes"]
samples = config["samples"]

rule all:
    input:
        expand("output/pmdv/{sample}/{sample}.vcf.gz",sample=samples),
        expand("output/nanopolish/{sample}/{sample}.methylation_frequency.tsv",sample=samples),
        expand("output/nanomonsv/{sample}/{sample}.nanomonsv.result.vcf",sample=samples),
        expand("output/whatshap/{sample}/{sample}.happlotagged_HP2.bam",sample=samples)

rule minimap2:
    input:
        fastq = lambda wildcards: config["fastq"][wildcards.sample],
    params:
        sh="script/minimap2/snakemake_minimap2.sh",
        ref=config["ref"],
        minimap_options=config["MINIMAP_OPTIONS"],
        samtools_options=config["SAMTOOLS_VIEW_OPTIONS"],
        samtools_sort=config["SAMTOOLS_SORT_OPTIONS"],
        samtools_index=config["SAMTOOLS_INDEX_OPTIONS"]
    output:
        "output/minimap2/{sample}/{sample}.bam"
    log:
        "log/minimap2/minimap2_{sample}_log"
    singularity:
        config["singularity_images"]["minimap2"]
    shell:
        "/bin/bash {params.sh} {params.ref} {input.fastq} {output} '{params.minimap_options}' '{params.samtools_options}' '{params.samtools_sort}' '{params.samtools_index}'"        

rule pmdv:
    input:
        "output/minimap2/{sample}/{sample}.bam"
    params:
        sh="script/pmdv/snakemake_pmdv.sh",
        ref=config["ref"],
        pmdv_threads=config["PMDV_THREADS"],
        pmdv_mode=config["PMDV_MODE"]
    output:
        "output/pmdv/{sample}/{sample}.vcf.gz"
    log:
        "log/pmdv/pmdv_{sample}_log"
    singularity:
        config["singularity_images"]["pmdv"]
    shell:
        "/bin/bash {params.sh} {input} {params.ref} {output} {params.pmdv_threads} {params.pmdv_mode}"

rule nanopolish:
    input:
        fast5_dir = lambda wildcards: config["fast5_dir"][wildcards.sample],
        fastq = lambda wildcards: config["fastq"][wildcards.sample],
        bam="output/minimap2/{sample}/{sample}.bam"
    params:
        sh="script/nanopolish/snakemake_nanopolish.sh",
        ref=config["ref"],
        hdf5_plugin_path=config["NANOPOLISH_HDF5_PLUGIN_PATH"],
        nanopolish_threads=config["NANOPOLISH_THREADS"]
    output:
        "output/nanopolish/{sample}/{sample}.methylation_frequency.tsv"
    log:
        "log/nanopolish_{sample}_log"
    singularity:
        config["singularity_images"]["nanopolish"]
    shell:
        "/bin/bash {params.sh} {params.ref} {params.hdf5_plugin_path} {input} {output} {params.nanopolish_threads}"

rule nanomonsv:
    input:
        "output/minimap2/{sample}/{sample}.bam"
    params:
        sh="script/nanomonsv/snakemake_nanomonsv.sh",
        ref=config["ref"],
        nanomonsv_options=config["NANOMONSV_OPTIONS"],
        control_panel=config["NANOMONSV_CONTROL_PANEL"],
        control_bam=config["NANOMONSV_CONTROL_BAM"],
        control_prefix=config["NANOMONSV_CONTROL_PREFIX"] 
    output:
        "output/nanomonsv/{sample}/{sample}.nanomonsv.result.vcf"
    log:
        "log/nanomonsv_{sample}_log"
    singularity:
        config["singularity_images"]["nanomonsv"]
    shell:
        "/bin/bash {params.sh} {params.ref} {input} {output} {params.control_panel} {params.control_bam} {params.control_prefix} '{params.nanomonsv_options}'"

rule glimpse:
    input:
        "output/minimap2/{sample}/{sample}.bam"
    params:
        sh="script/glimpse/snakemake_glimpse.sh",
        ref=config["ref"],
        panel_vcf_dir=config["GLIMPSE_PANEL_DIR"],
        map_file="/tools/GLIMPSE/maps/genetic_maps.b38/chr{chr_num}.b38.gmap.gz"
    output:
        "output/glimpse/{sample}/{sample}.chr{chr_num}.phased.bcf"
    log:
        "log/glimpse_{sample}_chr{chr_num}.log"
    singularity:
        config["singularity_images"]["glimpse"]
    shell:
        "/bin/bash {params.sh} {params.ref} {params.panel_vcf_dir} {input} {output} {params.map_file}"

rule glimpse_concat:
    input:
        lambda wildcards: expand("output/glimpse/{sample}/{sample}.chr{chr_num}.phased.bcf",sample=wildcards.sample,chr_num=chr_nums)
    params:
        sh="script/glimpse/snakemake_concat.sh"
    output:
        "output/glimpse/{sample}/{sample}.allchr.concat.vcf.gz"
    log:
        "log/glimpse_concat_{sample}_log"
    singularity:
        config["singularity_images"]["minimap2"]
    shell:
        "bash {params.sh} {output} {input}"

rule whatshap_haplotag:
    input:
        vcf="output/glimpse/{sample}/{sample}.allchr.concat.vcf.gz",
        bam="output/minimap2/{sample}/{sample}.bam"
    params:
        sh="script/whatshap/snakemake_whatshap_haplotag.sh",
        ref=config["ref"],
        whatshap_options=config["WHATSHAP_OPTIONS"],
        whatshap_threads=config["WHATSHAP_THREADS"]
    output:
        "output/whatshap/{sample}/{sample}.happlotagged.bam"
    log:
        "log/whatshap_haplotag_{sample}_log"
    singularity:
        config["singularity_images"]["whatshap"]
    shell:
        "/bin/bash {params.sh} {params.ref} {input.bam} {input.vcf} {output} '{params.whatshap_options}' '{params.whatshap_threads}'"

rule whatshap_split:
    input:
        "output/whatshap/{sample}/{sample}.happlotagged.bam"
    params:
        sh="script/whatshap/snakemake_whatshap_split.sh",
        whatshap_samtools_options=config["WHATSHAP_SAMTOOLS_OPTIONS"]
    output:
        HP1_bam="output/whatshap/{sample}/{sample}.happlotagged_HP1.bam",
        HP2_bam="output/whatshap/{sample}/{sample}.happlotagged_HP2.bam"
    log:
        "log/whatshap_split_{sample}_log"
    singularity:
        config["singularity_images"]["minimap2"]
    shell:
        "/bin/bash {params.sh} {input} {output.HP1_bam} {output.HP2_bam} '{params.whatshap_samtools_options}'"
