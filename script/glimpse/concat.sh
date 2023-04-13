#!/bin/sh

set -x
set -o errexit
set -o nounset
set -o pipefail

SAMPLE=$1
OUTPUT_DIR=$2

bcftools concat -Oz ${OUTPUT_DIR}/${SAMPLE}.chr*.phased.bcf > ${OUTPUT_DIR}/${SAMPLE}.allchr.concat.vcf.gz
tabix -p vcf -f ${OUTPUT_DIR}/${SAMPLE}.allchr.concat.vcf.gz
