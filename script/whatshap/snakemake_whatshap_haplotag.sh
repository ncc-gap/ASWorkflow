#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

REF=$1
INPUT_BAM=$2 #output/bam/genebay-sample-202107-17-03-1_CMK_3rd_mini.bam
INPUT_PHASED_VCF=$3 #"output/glimpse/genebay-sample-202107-17-03-1_CMK_3rd.allchr.concat.vcf.gz"
OUTPUT_BAM=$4
WHATSHAP_OPTIONS=$5
WHATSHAP_THREADS=$6

OUTPUT_DIR=`dirname ${OUTPUT_BAM}`

#haplotagged.BAM
whatshap haplotag \
${WHATSHAP_OPTIONS} \
--output-haplotag-list ${OUTPUT_DIR}/haplotag_list.tsv \
--output-threads ${WHATSHAP_THREADS} \
-o ${OUTPUT_BAM} \
--reference ${REF} \
${INPUT_PHASED_VCF} \
${INPUT_BAM}

