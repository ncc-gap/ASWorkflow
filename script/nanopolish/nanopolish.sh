#!/bin/bash -x

set -x
set -o errexit
set -o nounset

REF=$1
SAMPLE=$2
HDF5_PLUGIN_PATH=$3
FAST5_DIR=$4
INPUT_FASTQ=$5
INPUT_BAM=$6
OUTPUT_DIR=$7
NANOPOLISH_THREADS=$8

export HDF5_PLUGIN_PATH=${HDF5_PLUGIN_PATH}

#index
/tools/nanopolish/nanopolish index \
-d ${FAST5_DIR} ${INPUT_FASTQ}

#call-methylation
/tools/nanopolish/nanopolish call-methylation \
-t ${NANOPOLISH_THREADS} \
-r ${INPUT_FASTQ} \
-b ${INPUT_BAM} \
-g ${REF} \
> ${OUTPUT_DIR}/${SAMPLE}.methylation_calls.tsv

#convert2tsv
python ./script/nanopolish/calculate_methylation_frequency.py \
${OUTPUT_DIR}/${SAMPLE}.methylation_calls.tsv \
-s \
> ${OUTPUT_DIR}/${SAMPLE}.methylation_frequency.tsv
