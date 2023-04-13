#! /bin/bash

set -x
set -o errexit
set -o nounset
set -o pipefail

PMDV_IMAGE=$1
INPUT_BAM=$2
REF=$3
OUTPUT_DIR=$4
OUTPUT_PREFIX=$5
THREADS=$6
PMDV_MODE=$7

#PEPPER_MARGIN_DEEPVARIANT
singularity exec --bind /usr/lib/locale/ ${PMDV_IMAGE} \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_BAM}" \
-f "${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t "${THREADS}" \
${PMDV_MODE}

#pmdv.sh

