#!/bin/bash

set -x
set -o errexit
set -o nounset
set -o pipefail

REF=$1
INPUT_BAM=$2
OUTPUT_FILE=$3  
CONTROL_PANEL=$4
CONTROL_BAM=$5
CONTROL_PREFIX=$6
NANOMONSV_OPTIONS=$7

OUTPUT_PREFIX=`dirname ${OUTPUT_FILE}`/`basename ${OUTPUT_FILE} .nanomonsv.result.vcf` 

#nanomonsv_parse
nanomonsv parse ${INPUT_BAM} \
${OUTPUT_PREFIX}

#nanomonsv_get
nanomonsv get \
${OUTPUT_PREFIX} \
${INPUT_BAM} \
${REF} \
${NANOMONSV_OPTIONS} \
--control_panel_prefix ${CONTROL_PANEL} \
--control_prefix ${CONTROL_PREFIX} \
--control_bam ${CONTROL_BAM}

