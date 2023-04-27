#!/bin/sh

set -xv
set -o errexit
set -o nounset
set -o pipefail

FASTQ_PATH=$1
FAST5_DIR=$2
OUTPUT_DIR=$3
SOURCE_FILE=$4
SAMPLE=$5

source ${SOURCE_FILE}

#MINIMAP2
#required memory 20Gx, 30G◯
mkdir -p ${OUTPUT_DIR}/minimap2/${SAMPLE}

singularity exec \
${MINIMAP_IMAGE} \
/bin/bash -xv ./script/minimap2/minimap2.sh \
${REF} \
${SAMPLE} \
${FASTQ_PATH} \
${OUTPUT_DIR}/minimap2/${SAMPLE} \
"${MINIMAP_OPTIONS}" \
"${SAMTOOLS_VIEW_OPTIONS}" \
"${SAMTOOLS_SORT_OPTIONS}" \
"${SAMTOOLS_INDEX_OPTIONS}"

#PMDV
#required memory 20G_×:1thread、120G_◯(8threads)
mkdir -p ${OUTPUT_DIR}/pmdv/${SAMPLE}

/bin/bash -xv ./script/pmdv/pmdv.sh \
${PMDV_IMAGE} \
${OUTPUT_DIR}/minimap2/${SAMPLE}/${SAMPLE}_sorted.bam \
${REF} \
${OUTPUT_DIR}/pmdv/${SAMPLE} \
${SAMPLE} \
${PMDV_THREADS} \
${PMDV_MODE}

#NANOMONSV
mkdir -p ${OUTPUT_DIR}/nanomonsv/${SAMPLE}
mkdir -p ${OUTPUT_DIR}/nanomonsv/control

singularity exec \
${NANOMONSV_IMAGE} \
/bin/bash -xv ./script/nanomonsv/nanomonsv.sh \
${REF} \
${CONTROL_BAM} \
${OUTPUT_DIR}/nanomonsv/control/control \
${OUTPUT_DIR}/minimap2/${SAMPLE}/${SAMPLE}_sorted.bam \
${OUTPUT_DIR}/nanomonsv/${SAMPLE}/${SAMPLE} \
"${NANOMONSV_OPTIONS}" \
${CONTROL_PANEL}

#NANOPOLISH
#memory 10g
mkdir -p ${OUTPUT_DIR}/nanopolish/${SAMPLE}
singularity exec \
${NANOPOLISH_IMAGE} \
/bin/bash -xv ./script/nanopolish/nanopolish.sh \
${REF} \
${SAMPLE} \
${HDF5_PLUGIN_PATH} \
${FAST5_DIR} \
${INPUT_FASTQ} \
${OUTPUT_DIR}/minimap2/${SAMPLE}/${SAMPLE}_sorted.bam \
${OUTPUT_DIR}/nanopolish/${SAMPLE} \
${NANOPOLISH_THREADS}

#GLIMPSE
#memory 10g
mkdir -p ${OUTPUT_DIR}/glimpse/${SAMPLE}

cat ${CHR_NUM_TEXT} | while read line;
do
singularity exec \
${GLIMPSE_IMAGE} \
/bin/bash -xv ./script/glimpse/glimpse.sh \
${REF} \
${SAMPLE} \
${PANEL_VCF_DIR} \
${OUTPUT_DIR}/minimap2/${SAMPLE}/${SAMPLE}_sorted.bam \
${OUTPUT_DIR}/glimpse/${SAMPLE} \
${MAP_DIR} \
${line}
done

singularity exec \
${MINIMAP_IMAGE} \
/bin/bash -xv ./script/glimpse/concat.sh \
${SAMPLE} \
${OUTPUT_DIR}/glimpse/${SAMPLE}

#WHATSHAP
#required memory 10g
mkdir -p ${OUTPUT_DIR}/whatshap/${SAMPLE}

singularity exec \
${WHATSHAP_IMAGE} \
/bin/bash -xv ./script/whatshap/whatshap_haplotag.sh \
${REF} \
${SAMPLE} \
${OUTPUT_DIR}/minimap2/${SAMPLE}/${SAMPLE}_sorted.bam \
${OUTPUT_DIR}/glimpse/${SAMPLE}/${SAMPLE}.allchr.concat.vcf.gz \
${OUTPUT_DIR}/whatshap/${SAMPLE} \
"${WHATSHAP_OPTIONS}" \
${WHATSHAP_THREADS} 

singularity exec \
${MINIMAP_IMAGE} \
/bin/bash -xv ./script/whatshap/whatshap_split.sh \
${SAMPLE} \
${OUTPUT_DIR}/whatshap/${SAMPLE} \
"${WHATSHAP_SAMTOOLS_OPTIONS}"

