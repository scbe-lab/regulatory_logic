#!/bin/bash
BOWTIE1_REF="/omics-storage/lab-solana/genomes/Smed/Rink/schMedS3/bowtie1_colorspace_index/schMedS3_h1.fa.bowtie1"
OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/prep_kao/bowtie1C_output/"
SAMPLE_FOLDER="/omics-storage/lab-solana/alberto/DATA/static/reads/prep_kao/"

for f in ${SAMPLE_FOLDER}/* ; do

    if [ ! -d $f ] ; then
        continue
    fi

    x=${f##*/}
    z=${x}
    echo Starting with sample ${z} ...
    mkdir -p ${OUTPUT_FOLDER}/bowtie1C_out_${z}

    bowtie -C -S -p 12 \
        ${BOWTIE1_REF} \
        ${f}/*.fastq.gz | \
        samtools view -b - \
        > ${OUTPUT_FOLDER}/bowtie1C_out_${z}/${z}.bam

    echo "done with sample ${z} ..."

done
echo "Done."
