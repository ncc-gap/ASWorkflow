#!/bin/bash -x

set -x
set -o errexit
set -o nounset
set -o pipefail


REF=$1
HDF5_PLUGIN_PATH=$2
FAST5_DIR=$3
INPUT_FASTQ=$4
INPUT_BAM=$5
OUTPUT_TSV=$6 #"output/nanopolish/genebay-sample-202107-17-03-1_CMK_3rd_mini.methylation_frequency.tsv"
NANOPOLISH_THREADS=$7

OUTPUT_DIR=`dirname ${OUTPUT_TSV}`

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
> ${OUTPUT_DIR}/methylation_calls.tsv

#convert2tsv
python ./script/nanopolish/calculate_methylation_frequency.py \
${OUTPUT_DIR}/methylation_calls.tsv \
-s \
> ${OUTPUT_TSV}

