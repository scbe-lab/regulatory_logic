#!/bin/bash
KALLISTO_INDEX="/omics-storage/lab-solana/genomes/Smed/Rink/schMedS3/transdecoder_and_annotation/transdecoder/transcripts.fna.kallisto.index"
SAMPLES_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/data/bulk_rna/hnf4_combined_knockdown_bulk/"
KALLISTO_OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/hnf4_combined_knockdown_bulk/kallisto_output/raw_kallisto_output"

echo "output folder is $KALLISTO_OUTPUT_FOLDER"

# create output folder if it does not exist
if [ ! -d $f ] ; then
    mkdir -p $KALLISTO_OUTPUT_FOLDER
fi

# all the single pair-ends
for i in ${SAMPLES_FOLDER}/* ; do

    x=${i##*/}
    echo Starting with sample ${x} ...

    R1=${i}/*1.fq.gz
    R2=${i}/*2.fq.gz

    mkdir -p ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}

    kallisto quant -t 12 -i \
    $KALLISTO_INDEX \
    -o ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x} \
    $R1 $R2
    echo "done with sample ${x} ..."

done

echo "Done."
