#!/bin/sh

set -x
set -o errexit
set -o nounset
set -o pipefail

REF=$1
SAMPLE=$2
FASTQ_PATH=$3
MINIMAP_OUTPUT_DIR=$4
MINIMAP_OPTIONS=$5
SAMTOOLS_VIEW_OPTIONS=$6
SAMTOOLS_SORT_OPTIONS=$7
SAMTOOLS_INDEX_OPTIONS=$8


OUTPUT_BAM=${MINIMAP_OUTPUT_DIR}/${SAMPLE}_sorted.bam

minimap2 ${MINIMAP_OPTIONS} ${REF} ${FASTQ_PATH} | samtools view ${SAMTOOLS_VIEW_OPTIONS} > ${OUTPUT_BAM}.unsorted

samtools sort ${SAMTOOLS_SORT_OPTIONS} ${OUTPUT_BAM}.unsorted -o ${OUTPUT_BAM}

samtools index $SAMTOOLS_INDEX_OPTIONS ${OUTPUT_BAM}

rm ${OUTPUT_BAM}.unsorted

#minimap2.sh

