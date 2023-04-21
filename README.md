# ASWorkflow　(Adaptive_Sampling_Workflow)

## In preparation

## Introduction

Workflow for adaptive sampling

### How to setup

1.Preperation

git clone https://github.com/ncc-gap/ASWorkflow.git

cd ASWorkflow

2.Download files

`bash download_files.sh`

3.Pull singularity image

All singularity images are stored at ...

4.Make config file

Edit reference files, image files, and more.

`vi ${config}`

### How to run
`bash run.sh ${input_fastq_path} ${input_fast5_dir} ${output_dir} ${config} ${sample_id}`



## Snakemake
### How to setup

Preparation

### Pull singularity images
```
mkdir $PWD/image
singularity pull $PWD/image/pepper_deepvariant_r0.8.sif docker://kishwars/pepper_deepvariant:r0.8
singularity pull $PWD/image/nanomonsv_v0.5.0.sif docker://friend1ws/nanomonsv:v0.5.0
singularity pull $PWD/image/glimpse_v1.1.1.sif docker://wn505/glimpse:v1.1.1
singularity pull $PWD/image/minimap2_v2.22_2.sif docker://wn505/minimap2:v2.22_2
singularity pull $PWD/image/whatshap_v1.4.0.sif docker://wn505/whatshap:v1.4.0
singularity pull $PWD/image/nanopolish_v0.13.3.sif docker://aokad/nanopolish:0.0.1
```

### nanopolish
### Download hd5_plugin
```
mkdir -p $PWD/downloads/nanopolish
wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz \
-O $PWD/downloads/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
tar -zxvf $PWD/downloads/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz -C $PWD/downloads/nanopolish/
HDF5_PLUGIN_PATH=$PWD/downloads/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin
```

### nanomonsv
### Download control_panel
```
mkdir -p $PWD/downloads/nanomonsv/control_panel
wget https://zenodo.org/api/files/5c116b75-6ef0-4445-9fa8-c5989639da5f/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz \
-O $PWD/downloads/nanomonsv/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz
tar -xvf $PWD/downloads/nanomonsv/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz -C $PWD/downloads/nanomonsv/control_panel/
NANOMONSV_CONTROL_PANEL=$PWD/downloads/nanomonsv/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control
```

### Download fastq
```
mkdir -p $PWD/downloads/nanomonsv/control_bam
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20211013_ONT_Rebasecalled/NA18989/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom_fastq_pass.fastq.gz \
-O $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom_fastq_pass.fastq.gz
```

### Download reference
```
mkdir -p $PWD/reference
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
  -O $PWD/reference/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai \
  -O $PWD/reference/Homo_sapiens_assembly38.fasta.fai
```

### Alignment(Make bam file)
#### Takes 24 hours even with qsub
```
singularity exec $PWD/image/minimap2_v2.22_2.sif sh -c \
  "minimap2 -ax map-ont -t 8 -p 0.1 $PWD/reference/Homo_sapiens_assembly38.fasta $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom_fastq_pass.fastq.gz \
  | samtools view -Shb > $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.unsorted.bam && \
  samtools sort -@ 8 -m 2G $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.unsorted.bam \
  -o $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.bam && \
  samtools index $PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.bam"
NANOMONSV_CONTROL_BAM=$PWD/downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.bam
```

### Make control prefix
```
mkdir -p $PWD/downloads/nanomonsv/control_prefix
singularity exec $PWD/image/nanomonsv_v0.5.0.sif sh -c \
  "nanomonsv parse downloads/nanomonsv/control_bam/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom.bam \
   $PWD/downloads/nanomonsv/control_prefix/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom"
NANOMONSV_CONTROL_PREFIX=$PWD/downloads/nanomonsv/control_prefix/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom
```

### Glimpse
### Download conrtol_panel and make vcf
```
mkdir -p $PWD/downloads/glimpse 
for i in {1..22} 
do
    PANEL_VCF=$PWD/downloads/glimpse/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz
    REFBCF=$PWD/downloads/glimpse/1000GP.chr${i}.noNA12878.bcf
    REFVCF=$PWD/downloads/glimpse/1000GP.chr${i}.noNA12878.vcf.gz
    REFTSV=$PWD/downloads/glimpse/1000GP.chr${i}.noNA12878.tsv.gz
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz{,.tbi} \
    -O  ${PANEL_VCF}
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
