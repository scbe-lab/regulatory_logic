#!/bin/bash
ANNOT="/omics-storage/lab-solana/standard_references/Smed_Rink/Smed_Rink_MASKED.gtf"
OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/prep_kao/bowtie1C_htseq_counts/"
SAMPLE_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/prep_kao/bowtie1C_output/"

for f in ${SAMPLE_FOLDER}/* ; do

    if [ ! -d $f ] ; then
        continue
    fi

    x=${f##*/}
    z=${x//bowtie1C_out_}
    echo Starting with sample ${z} ...
    mkdir -p ${OUTPUT_FOLDER}/htseq_out_${z}

    htseq-count --nonunique all -f bam \
    -r pos --type gene \
    -s no -n 12 \
    ${f}/*.bam \
    $ANNOT \
    > ${OUTPUT_FOLDER}/htseq_out_${z}/${z}.counts \


    echo "done with sample ${z} ..."

done
echo "Done."
