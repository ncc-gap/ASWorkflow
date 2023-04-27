#!/bin/bash

set -x
set -o errexit
set -o nounset
set -o pipefail

input_control_bam=$1 

#4.1.1 Download hd5_plugin
mkdir -p $PWD/downloads/nanopolish
wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz \
-O $PWD/downloads/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
tar -zxvf $PWD/downloads/nanopolish/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz -C $PWD/downloads/nanopolish/

#4.2.1 Download reference
mkdir -p $PWD/reference
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
  -O $PWD/reference/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai \
  -O $PWD/reference/Homo_sapiens_assembly38.fasta.fai

#4.2.2 Download control_panel
mkdir -p $PWD/downloads/nanomonsv/control_panel
wget https://zenodo.org/api/files/5c116b75-6ef0-4445-9fa8-c5989639da5f/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz \
-O $PWD/downloads/nanomonsv/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz
tar -xvf $PWD/downloads/nanomonsv/control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control.tar.gz -C $PWD/downloads/nanomonsv/control_panel/

#4.2.5 Make control prefix
mkdir -p $PWD/downloads/nanomonsv/control_prefix
singularity exec $PWD/image/nanomonsv_v0.5.0.sif sh -c \
  "nanomonsv parse ${input_control_bam} \
   $PWD/downloads/nanomonsv/control_prefix/TP53_range_5M_THP.sorted"

#4.3.1 Download conrtol_panel and make vcf
mkdir -p $PWD/downloads/glimpse 
    PANEL_VCF=$PWD/downloads/glimpse/CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz
    REFBCF=$PWD/downloads/glimpse/1000GP.chr17.noNA12878.bcf
    REFVCF=$PWD/downloads/glimpse/1000GP.chr17.noNA12878.vcf.gz
    REFTSV=$PWD/downloads/glimpse/1000GP.chr17.noNA12878.tsv.gz
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz \
    -O  ${PANEL_VCF}
    singularity exec ../image/minimap2_v2.22_2.sif sh -c \
    "tabix -p vcf ${PANEL_VCF}"
    singularity exec ../image/glimpse_v1.1.1.sif sh -c \
    "bcftools norm -m -any ${PANEL_VCF} -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps -s ^NA12878,NA12891,NA12892 --threads 4 -Ob -o ${REFBCF} && \
    bcftools index -f ${REFBCF} --threads 4 && \
    bcftools norm -m -any ${PANEL_VCF} -Ou --threads 4 | bcftools view -G -m 2 -M 2 -v snps -Oz -o ${REFVCF} && \
    bcftools index -f ${REFVCF} --threads 4 && \
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${REFVCF} | bgzip -c > ${REFTSV} && \
    tabix -s1 -b2 -e2 ${REFTSV}"
