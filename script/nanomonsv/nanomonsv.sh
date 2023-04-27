#!/bin/sh

REF=$1
CONTROL_BAM=$2
CONTROL_PREFIX=$3
INPUT_BAM=$4
OUTPUT_DIR=$5
NANOMONSV_OPTIONS=$6
CONTROL_PANEL=$7

#control_parse
nanomonsv parse ${CONTROL_BAM} \
${CONTROL_PREFIX}

#nanomonsv_parse
nanomonsv parse ${INPUT_BAM} \
${OUTPUT_DIR}

#nanomonsv_get
nanomonsv get \
${OUTPUT_DIR} \
${INPUT_BAM} \
${REF} \
${NANOMONSV_OPTIONS} \
--control_panel_prefix ${CONTROL_PANEL} \
--control_prefix ${CONTROL_PREFIX} \
--control_bam ${CONTROL_BAM}

