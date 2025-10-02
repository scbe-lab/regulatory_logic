#!/bin/bash
BOWTIE2_INDEX="/omics-storage/lab-solana/genomes/Smed/Rink/schMedS3/schMedS3_h1.bowtie2"
BOWTIE2_OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/ap_stuckemann/bowtie2_output/"
SAMPLE_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/ap_stuckemann/cutadapt_output"

for f in ${SAMPLE_FOLDER}/* ; do

    if [ ! -d $f ] ; then
        continue
    fi

    x=${f##*/}
    z=${x//cutadapt_out_}
    echo Starting with sample ${z} ...
    mkdir -p ${BOWTIE2_OUTPUT_FOLDER}/bowtie2_out_${z}

    bowtie2 -x $BOWTIE2_INDEX \
    -r $f/*.gz \
    --quiet -p 20 \
    -X 75 | \
    samtools view -b - | samtools sort -m 10G -@ 12 > ${BOWTIE2_OUTPUT_FOLDER}/bowtie2_out_${z}/${z}.bam
    echo "done with sample ${z} ..."

done
echo "Done."
