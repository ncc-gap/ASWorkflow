#!/bin/bash

set -x
set -o errexit
set -o nounset
set -o pipefail

REF=$1 
PANEL_DIR=$2
INPUT_BAM=$3
PHASED_BCF=$4
MAP=$5

OUTPUT_DIR=`dirname ${PHASED_BCF}`
CONTROL_DIR=`dirname ${OUTPUT_DIR}`/control

SAMPLE=`basename ${OUTPUT_DIR}`

CHR_NUM_WITH_TAIL=${PHASED_BCF##*chr}
CHR_NUM=${CHR_NUM_WITH_TAIL%%.phased.bcf}

PANEL_VCF=${PANEL_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR_NUM}.filtered.shapeit2-duohmm-phased.vcf.gz

REFBCF=${PANEL_DIR}/1000GP.chr${CHR_NUM}.noNA12878.bcf
REFVCF=${PANEL_DIR}/1000GP.chr${CHR_NUM}.noNA12878.vcf.gz
REFTSV=${PANEL_DIR}/1000GP.chr${CHR_NUM}.noNA12878.tsv.gz

INTER_MPILEUP_VCF=${OUTPUT_DIR}/${SAMPLE}_mpileup.chr${CHR_NUM}.vcf.gz

#if [ ! -f ${REFBCF} ] ;then
#bcftools norm -m -any ${PANEL_VCF} -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps -s ^NA12878,NA12891,NA12892 --threads 4 -Ob -o ${REFBCF}
#bcftools index -f ${REFBCF} --threads 4
#fi

#if [ ! -f ${REFVCF} ] ;then
#bcftools norm -m -any ${PANEL_VCF} -Ou --threads 4 | bcftools view -G -m 2 -M 2 -v snps -Oz -o ${REFVCF}
#bcftools index -f ${REFVCF} --threads 4
#fi

#if [ ! -f ${REFTSV} ] ;then
#bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${REFVCF} | bgzip -c > ${REFTSV}
#tabix -s1 -b2 -e2 ${REFTSV}
#fi

bcftools mpileup -f ${REF} -X ont -a 'FORMAT/DP' -T ${REFVCF} -r chr${CHR_NUM} ${INPUT_BAM} -Ou | bcftools call -P 0.01 -Aim -C alleles -T ${REFTSV} -Oz -o ${INTER_MPILEUP_VCF}
bcftools index -f ${INTER_MPILEUP_VCF}

#4.1 glimpse
GLIMPSE_chunk --input ${REFVCF} --region chr${CHR_NUM} --window-size 2000000 --buffer-size 200000 --output ${OUTPUT_DIR}/1000GP.chr${CHR_NUM}.chunks.txt

#5.1 glimpse
while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=${OUTPUT_DIR}/${SAMPLE}.chr${CHR_NUM}.${ID}.bcf
GLIMPSE_phase  --thread 8 --input ${INTER_MPILEUP_VCF} --reference ${REFBCF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
bcftools index -f ${OUT}
done < ${OUTPUT_DIR}/1000GP.chr${CHR_NUM}.chunks.txt

#6.1 glimpse_ligate
LST=${OUTPUT_DIR}/${SAMPLE}.list.chr${CHR_NUM}.txt
ls ${OUTPUT_DIR}/${SAMPLE}.chr${CHR_NUM}.*.bcf > ${LST}
MERGED_BCF=${OUTPUT_DIR}/${SAMPLE}.chr${CHR_NUM}.merged.bcf
GLIMPSE_ligate --input ${LST} --output ${MERGED_BCF}
bcftools index -f ${MERGED_BCF}

#7.1 glimpse_sample
GLIMPSE_sample --input ${MERGED_BCF} --solve --output ${PHASED_BCF}
bcftools index -f ${PHASED_BCF}

rm ${OUTPUT_DIR}/${SAMPLE}.chr${CHR_NUM}.[0-9][0-9].bcf
rm ${OUTPUT_DIR}/${SAMPLE}.chr${CHR_NUM}.[0-9][0-9].bcf.csi
rm ${OUTPUT_DIR}/${SAMPLE}.list.chr${CHR_NUM}.txt
rm ${MERGED_BCF}
rm ${MERGED_BCF}.csi
rm ${INTER_MPILEUP_VCF}
rm ${INTER_MPILEUP_VCF}.csi
rm ${OUTPUT_DIR}/1000GP.chr${CHR_NUM}.chunks.txt
