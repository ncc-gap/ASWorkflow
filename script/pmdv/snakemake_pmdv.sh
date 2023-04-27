#! /bin/bash

set -x
set -o errexit
set -o nounset
set -o pipefail

INPUT_BAM=$1
REF=$2
OUTPUT_VCF=$3
THREADS=$4
PMDV_MODE=$5

OUTPUT_DIR=`dirname ${OUTPUT_VCF}`
OUTPUT_PREFIX=`basename ${OUTPUT_VCF} .vcf.gz`

run_pepper_margin_deepvariant call_variant \
-b "${INPUT_BAM}" \
-f "${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t "${THREADS}" \
${PMDV_MODE}
