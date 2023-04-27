#!/bin/sh

set -x
set -o errexit
set -o nounset
set -o pipefail

REF=$1
SAMPLE=$2
PANEL_VCF_DIR=$3
INPUT_BAM=$4
OUTPUT_DIR=$5
MAP_DIR=$6
CHR_NUM=$7

PANEL_VCF=${PANEL_VCF_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR_NUM}.filtered.shapeit2-duohmm-phased_filter_HWE_0.000001.vcf.gz
REFBCF=${OUTPUT_DIR}/1000GP.chr${CHR_NUM}.noNA12878.bcf
REFVCF=${OUTPUT_DIR}/1000GP.chr${CHR_NUM}.noNA12878.vcf.gz
REFTSV=${OUTPUT_DIR}/1000GP.chr${CHR_NUM}.noNA12878.tsv.gz
INTER_MPILEUP_VCF=${OUTPUT_DIR}/${SAMPLE}_mpileup.chr${CHR_NUM}.vcf.gz

bcftools norm -m -any ${PANEL_VCF} -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps -s ^NA12878,NA12891,NA12892 --threads 4 -Ob -o ${REFBCF}
bcftools index -f ${REFBCF} --threads 4

bcftools norm -m -any ${PANEL_VCF} -Ou --threads 4 | bcftools view -G -m 2 -M 2 -v snps -Oz -o ${REFVCF}
bcftools index -f ${REFVCF} --threads 4

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${REFVCF} | bgzip -c > ${REFTSV}
tabix -s1 -b2 -e2 ${REFTSV}

bcftools mpileup -f ${REF} -X ont -a 'FORMAT/DP' -T ${REFVCF} -r chr${CHR_NUM} ${INPUT_BAM} -Ou | bcftools call -P 0.01 -Aim -C alleles -T ${REFTSV} -Oz -o ${INTER_MPILEUP_VCF}
bcftools index -f ${INTER_MPILEUP_VCF}

#4.1 glimpse
GLIMPSE_chunk --input ${REFVCF} --region chr${CHR_NUM} --window-size 2000000 --buffer-size 200000 --output ${OUTPUT_DIR}/1000GP.chr${CHR_NUM}.chunks.txt

#5.1 glimpse
MAP=${MAP_DIR}/chr${CHR_NUM}.b38.gmap.gz
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
PHASED_BCF=${OUTPUT_DIR}/${SAMPLE}.chr${CHR_NUM}.phased.bcf
GLIMPSE_sample --input ${MERGED_BCF} --solve --output ${PHASED_BCF}
bcftools index -f ${PHASED_BCF}
