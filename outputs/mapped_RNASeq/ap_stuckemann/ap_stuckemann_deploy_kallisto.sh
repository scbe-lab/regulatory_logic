#!/bin/bash
KALLISTO_INDEX="/omics-storage/lab-solana/genomes/Smed/Rink/schMedS3/transdecoder_and_annotation/transdecoder/transcripts.fna.kallisto.index"
KALLISTO_OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/ap_stuckemann/kallisto_output/raw_kallisto_output"
SAMPLE_FOLDER="/omics-storage/lab-solana/alberto/DATA/static/reads/smed_ap_stuckemann"

for f in ${SAMPLE_FOLDER}/* ; do

    if [ ! -d $f ] ; then
        continue
    fi

    x=${f##*/}
    echo Starting with sample ${x} ...
    mkdir -p ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}
    kallisto quant -t 12 -i \
    $KALLISTO_INDEX \
    --single $f/*.f*.gz \
    -l 76 -s 10 \
    -o ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}
    echo "done with sample ${x} ..."

done
echo "Done."
