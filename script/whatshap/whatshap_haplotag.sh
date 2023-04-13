#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

REF=$1
SAMPLE=$2
INPUT_BAM=$3
INPUT_PHASED_BCF=$4
OUTPUT_DIR=$5
WHATSHAP_OPTIONS=$6
WHATSHAP_THREADS=$7

#haplotagged.BAM
whatshap haplotag \
${WHATSHAP_OPTIONS} \
--output-haplotag-list ${OUTPUT_DIR}/haplotag_list.tsv \
--output-threads ${WHATSHAP_THREADS} \
-o ${OUTPUT_DIR}/${SAMPLE}.happlotagged.bam \
--reference ${REF} \
${INPUT_PHASED_BCF} \
${INPUT_BAM}
