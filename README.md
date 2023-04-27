# ASWorkflow　(Adaptive_Sampling_Workflow)

## Introduction

This workflow is an efficient computational workflow for Adaptive Sampling.

This workflow consists of the following steps 
#### Note that Basecalling is assumed to be done separately. Please specify the paths to the FASTQ and FAST5 files individually in the config.

### 1. Read mapping
Map the reads to the reference genome using __Minimap2__.

### 2. Variant calling
Call the variants using __Pepper-Margin-DeepVariant__.

### 3. Structual Variation calling
Call the structual variation using __nanomonsv__.

### 4. Genotype imputation
Referring to the reference panel, low-coverage reads in off-target regions will also be used to perform genome-wide phasing by __GLIMPSE__.

### 5. Haplotype phasing
Tagging reads by haplotype for visualization using __Whatshap__.

### 6. Methylation calling
Determine methylation status using __Nanopolish__.


## Tutorial

This tutorial provides a step-by-step guide to analyze the sequence data (FASTQ and FAST5 files) of the CMK-86 cell line. 

The workflow is designed to cover the region surrounding TP53 gene.


### How to setup

__Comment__

Each analysis uses a __singularity image__, so please make sure you can use __Singularity__ beforehand.

__1. Download this repogitry__

we use the `convert_bam_for_methylation.py` script (https://github.com/timplab/nanopore-methylation-utilities) in ./script/nanopolish.

```
git clone https://github.com/ncc-gap/ASWorkflow.git
cd tutorial
```

__2. Install snakemake__

We also use __Snakemake__ to control the workflow, so please install __Snakemake__ beforehand.

Please see snakemake official document.

```
virtualenv -p python3 .venv
$ source .venv/bin/activate
$ pip install snakemake
```
or
```
pip3 install snakemake
```

__3. Pull singularity images__

__note that PWD is tutorial__

```
mkdir $PWD/image
singularity pull $PWD/image/pepper_deepvariant_r0.8.sif docker://kishwars/pepper_deepvariant:r0.8
singularity pull $PWD/image/nanomonsv_v0.5.0.sif docker://friend1ws/nanomonsv:v0.5.0
singularity pull $PWD/image/glimpse_v1.1.1.sif docker://wn505/glimpse:v1.1.1
singularity pull $PWD/image/minimap2_v2.22_2.sif docker://wn505/minimap2:v2.22_2
singularity pull $PWD/image/whatshap_v1.4.0.sif docker://wn505/whatshap:v1.4.0
singularity pull $PWD/image/nanopolish_v0.13.3.sif docker://aokad/nanopolish:0.0.1
```

__4. Download `input FAST5`, `input FASTQ` files and `control bam (for nanomonsv)` from zenodo__

```
Zenodo page is now preparing.
```

__5. Download tools and control data.__


__To obtain the tools and panel data needed to run the tutorial, you can use the following script__

Before executing the following commands, __please make sure that singularity is available__.

```
bash preparation.sh ${control_bam}
```

__This command will generate several files under the downloads directory.__

__Note that the following script uses the downloaded control data, so be sure to execute No.3(Pull singularity images) and No.4(Download `input FAST5`, `input FASTQ` files and `control bam (for nanomonsv)` from zenodo) before executing this No.5(Download tools and control data).__


__6.Additions to config__

The contents of config.cfg must be rewritten in __three places__.
1. input fastq file
2. fast5 directory of input
3. input control bam file (for nanomonsv)

```
vi config.cfg

1. Path of downloaded input_fastq

2. Path of downloaded input_fast5 (directory path)

3. Path of downloaded control_bam (for nanomonsv)
```


## For full-scale use, Please refer to the following "Download Tool, control_panels" to download further necessary data__

__Please see the `Download tools,control_panels` section at the bottom of this page__


### Final directory structure

```
.
├── REAME.md
├── Snakefile
├── config.yaml
├── preparetion.sh
├── downloads
│   ├── glimpse
│   │   ├──CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz
│   │   └──CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
│   ├── nanomonsv
│   │   ├── control_bam
│   │   │   ├── THP-1_control_tutorial_TP53_range_5M.sorted.bam (This file can be downloaded from zenodo and placed anywhere you like. 
│   │   │   │                                                                                       Please rewrite the placed location in your config.)
│   │   │   └── THP-1_control_tutorial_TP53_range_5M.sorted.bam.bai
│   │   ├── control_panel
│   │   │   ├── hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control
│   │   │   ├── ..........
│   │   │   └── hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.rearrangement.sorted.bedpe.gz.tbi
│   │   └── control_prefix
│   │   │   ├── THP-1_control_tutorial_TP53_range_5M.sorted.bed.gz
│   │   │   ├── ..........
│   │       └── THP-1_control_tutorial_TP53_range_5M.sorted.bedpe.gz.tbi
│   └── nanopolish
│       └──ont-vbz-hdf-plugin-1.0.1-Linux
├── image
│   ├── glimpse_v1.1.1.sif
│   ├── ..........
│   └── whatshap_v1.4.0.sif
├── input
│   ├── fast5_dir
│   │   └── CMK-86_input_tutorial_TP53_range_5M.fast5 (This file can be downloaded from zenodo and placed anywhere you like. 
│   │   │   │                                                                                       Please rewrite the placed location in your config.)
│   └── fastq
│       └── CMK-86_input_tutorial_TP53_range_5M.fastq (This file can be downloaded from zenodo and placed anywhere you like. 
│                                                                                        Please rewrite the placed location(directory) in your config.)
├── output
│   ├── glimpse
│   ├── minimap2
│   ├── nanomonsv
│   ├── nanopolish
│   ├── pmdv
│   └── whatshap
├── reference
│   ├── Homo_sapiens_assembly38.fasta
│   └── Homo_sapiens_assembly38.fasta.fai
└── script
    ├── glimpse
    │   ├── snakemake_concat.sh
    │   ├── snakemake_glimpse.sh
    ├── minimap2
    │   └── snakemake_minimap2.sh
    ├── nanomonsv
    │   └── snakemake_nanomonsv.sh
    ├── nanopolish
    │   ├── calculate_methylation_frequency.py
    │   ├── snakemake_nanopolish.sh
    ├── pmdv
    │   └── snakemake_pmdv.sh
    └── whatshap
        ├── snakemake_whatshap_haplotag.sh
<<<<<<< HEAD
        └──  snakemake_whatshap_split.sh

```


=======
        └── snakemake_whatshap_split.sh

```

>>>>>>> origin/main
## How to run

### For dry run
__Please not to forget to activate pyenv, if you have installed snakemake within pyenv.__

```
source {path/to/dir/bin/activate}

snakemake -np --verbose
```

if you create dag.png

```
snakemake --dag | dot -Tpng > dag.png
```

### For run
```
snakemake  --cores all --verbose --use-singularity
```

### Precise instruction of downloaing tools,control_panels.

__If you have already perfomed tutorial, you just need to execute 4.2.3(Download fastq), 4.2.5(Make control prefix) and 4.3.1(Download conrtol_panel and make vcf).__


#### __4.1 nanopolish__

__4.1.1 Download hd5_plugin__
```
mkdir -p $PWD/downloads/nanopolish
wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz \
-O $PWD/downloads/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
tar -zxvf $PWD/downloads/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz -C $PWD/downloads/nanopolish/
HDF5_PLUGIN_PATH=$PWD/downloads/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin
```

#### __4.2 nanomonsv__

__4.2.1 Download reference__
```
mkdir -p $PWD/reference
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
  -O $PWD/reference/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai \
  -O $PWD/reference/Homo_sapiens_assembly38.fasta.fai
```

__4.2.2 Download control_panel__
```
mkdir -p $PWD/downloads/nanomonsv/control_panel
wget https://zenodo.org/api/files/5c116b75-6ef0-4445-9fa8-c5989639da5f/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz \
-O $PWD/downloads/nanomonsv/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz
tar -xvf $PWD/downloads/nanomonsv/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz -C $PWD/downloads/nanomonsv/control_panel/
NANOMONSV_CONTROL_PANEL=$PWD/downloads/nanomonsv/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control
```

__4.2.3 Download fastq__
```
mkdir -p $PWD/downloads/nanomonsv/control_bam
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20211013_ONT_Rebasecalled/NA18989/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom_fastq_pass.fastq.gz \
-O $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom_fastq_pass.fastq.gz
```

__4.2.4 Alignment : Takes 24 hours even with qsub__
```
singularity exec $PWD/image/minimap2_v2.22_2.sif sh -c \
  "minimap2 -ax map-ont -t 8 -p 0.1 $PWD/reference/Homo_sapiens_assembly38.fasta $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom_fastq_pass.fastq.gz \
  | samtools view -Shb > $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.unsorted.bam && \
  samtools sort -@ 8 -m 2G $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.unsorted.bam \
  -o $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.bam && \
  samtools index $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.bam"
NANOMONSV_CONTROL_BAM=$PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.bam
```

__4.2.5 Make control prefix__
```
mkdir -p $PWD/downloads/nanomonsv/control_prefix
singularity exec $PWD/image/nanomonsv_v0.5.0.sif sh -c \
  "nanomonsv parse downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.bam \
   $PWD/downloads/nanomonsv/control_prefix/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom"
NANOMONSV_CONTROL_PREFIX=$PWD/downloads/nanomonsv/control_prefix/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom
```

#### __4.3 GLIMPSE__

__4.3.1 Download conrtol_panel and make vcf__
This script is based on the GLIMPSE1 tutorial but with slight modifications.

```
mkdir -p $PWD/downloads/glimpse 
for i in {1..22} 
do
    PANEL_VCF=$PWD/downloads/glimpse/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz
    REFBCF=$PWD/downloads/glimpse/1000GP.chr${i}.noNA12878.bcf
    REFVCF=$PWD/downloads/glimpse/1000GP.chr${i}.noNA12878.vcf.gz
    REFTSV=$PWD/downloads/glimpse/1000GP.chr${i}.noNA12878.tsv.gz
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz \
    -O  ${PANEL_VCF}
    singularity exec $PWD/image/minimap2_v2.22_2.sif sh -c \
    "tabix -p vcf ${PANEL_VCF}"
    singularity exec $PWD/image/glimpse_v1.1.1.sif sh -c \
    "bcftools norm -m -any ${PANEL_VCF} -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps -s ^NA12878,NA12891,NA12892 --threads 4 -Ob -o ${REFBCF} && \
    bcftools index -f ${REFBCF} --threads 4 && \
    bcftools norm -m -any ${PANEL_VCF} -Ou --threads 4 | bcftools view -G -m 2 -M 2 -v snps -Oz -o ${REFVCF} && \
    bcftools index -f ${REFVCF} --threads 4 && \
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${REFVCF} | bgzip -c > ${REFTSV} && \
    tabix -s1 -b2 -e2 ${REFTSV}"
done
GLIMPSE_PANEL_VCF_DIR=$PWD/downloads/glimpse
```
<<<<<<< HEAD
## License

ASWorkflow is free for academic use only. If you are not a member of a public funded academic and/or education and/or research institution you must obtain a commercial license from National Cancer Center; please email Yuichi Shiraishi (yuishira@ncc.go.jp). 

=======


## bash

1.Download this reoositry

```
git clone https://github.com/ncc-gap/ASWorkflow.git
cd ASWorkflow
```
2.Download files

`bash download_files.sh`

3.Pull singularity image

All singularity images are stored at ...

4.Make config file

Edit reference files, image files, and more.

`vi ${config}`
## How to run
`bash run.sh ${input_fastq_path} ${input_fast5_dir} ${output_dir} ${config} ${sample_id}`


## License

ASWorkflow is free for academic use only. If you are not a member of a public funded academic and/or education and/or research institution you must obtain a commercial license from National Cancer Center; please email Yuichi Shiraishi (yuishira@ncc.go.jp). 
>>>>>>> origin/main
