#!/bin/bash

set -x
set -o errexit
set -o nounset
set -o pipefail

OUTPUT_VCF=$1
INPUT_VCFS=${@:2}

bcftools concat -Oz ${INPUT_VCFS} > ${OUTPUT_VCF}

tabix -p vcf -f ${OUTPUT_VCF}
