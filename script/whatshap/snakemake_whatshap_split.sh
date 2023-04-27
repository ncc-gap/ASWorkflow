#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

INPUT_BAM=$1
OUTPUT_HP1_BAM=$2
OUTPUT_HP2_BAM=$3
WHATSHAP_SAMTOOLS_OPTIONS=$4

##index
samtools index \
${WHATSHAP_SAMTOOLS_OPTIONS} \
${INPUT_BAM}

##HP1
samtools view \
${INPUT_BAM} \
-d HP:1 \
-O bam \
-o ${OUTPUT_HP1_BAM} \
${WHATSHAP_SAMTOOLS_OPTIONS}

samtools index \
${WHATSHAP_SAMTOOLS_OPTIONS} \
${OUTPUT_HP1_BAM}

##HP2
samtools view \
${INPUT_BAM} \
-d HP:2 \
-O bam \
-o ${OUTPUT_HP2_BAM} \
${WHATSHAP_SAMTOOLS_OPTIONS}

samtools index \
${WHATSHAP_SAMTOOLS_OPTIONS} \
${OUTPUT_HP2_BAM}
