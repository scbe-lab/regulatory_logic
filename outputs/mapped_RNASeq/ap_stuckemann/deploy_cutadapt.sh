#!/bin/bash
CUTADAPT_OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/ap_stuckemann/cutadapt_output/"
SAMPLE_FOLDER="/omics-storage/lab-solana/alberto/DATA/static/reads/smed_ap_stuckemann"

for f in ${SAMPLE_FOLDER}/* ; do

    if [ ! -d $f ] ; then
        continue
    fi

    x=${f##*/}
    echo Starting with sample ${x} ...
    mkdir -p ${CUTADAPT_OUTPUT_FOLDER}/cutadapt_out_${x}
    cutadapt -a TCGTATGCCGTCT -n 5 -m 20 -q 25 -o ${CUTADAPT_OUTPUT_FOLDER}/cutadapt_out_${x}/${x}.trimmed.fq.gz -j 12 ${f}/*.fastq.gz
    echo "done with sample ${x} ..."

done
echo "Done."
