#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
OUTPUT_DIR=$2
WHATSHAP_SAMTOOLS_OPTIONS=$3

#index
samtools index \
${WHATSHAP_SAMTOOLS_OPTIONS} \
${OUTPUT_DIR}/${SAMPLE}.happlotagged.bam

#HP1
samtools view \
${OUTPUT_DIR}/${SAMPLE}.happlotagged.bam \
-d HP:1 \
-O bam \
-o ${OUTPUT_DIR}/${SAMPLE}.happlotagged_HP1.bam \
${WHATSHAP_SAMTOOLS_OPTIONS} 

samtools index \
${WHATSHAP_SAMTOOLS_OPTIONS} \
${OUTPUT_DIR}/${SAMPLE}.happlotagged_HP1.bam \

#HP2
samtools view \
${OUTPUT_DIR}/${SAMPLE}.happlotagged.bam \
-d HP:2 \
-O bam \
-o ${OUTPUT_DIR}/${SAMPLE}.happlotagged_HP2.bam \
${WHATSHAP_SAMTOOLS_OPTIONS} 

samtools index \
${WHATSHAP_SAMTOOLS_OPTIONS} \
${OUTPUT_DIR}/${SAMPLE}.happlotagged_HP2.bam

