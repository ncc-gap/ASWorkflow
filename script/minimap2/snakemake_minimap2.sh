#!/bin/bash

set -x
set -o errexit
set -o nounset
set -o pipefail

REF=$1
FASTQ_PATH=$2
OUTPUT_BAM=$3
MINIMAP_OPTIONS=$4
SAMTOOLS_VIEW_OPTIONS=$5
SAMTOOLS_SORT_OPTIONS=$6
SAMTOOLS_INDEX_OPTIONS=$7


minimap2 ${MINIMAP_OPTIONS} ${REF} ${FASTQ_PATH} | samtools view ${SAMTOOLS_VIEW_OPTIONS} > ${OUTPUT_BAM}.unsorted

samtools sort ${SAMTOOLS_SORT_OPTIONS} ${OUTPUT_BAM}.unsorted -o ${OUTPUT_BAM}

samtools index $SAMTOOLS_INDEX_OPTIONS ${OUTPUT_BAM}

rm ${OUTPUT_BAM}.unsorted
