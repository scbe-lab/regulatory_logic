#!/bin/bash
ANNOT="/omics-storage/lab-solana/standard_references/Smed_Rink/Smed_Rink_MASKED.gtf"

echo "NOW HTSEQ"

source /home/alberto/.bashrc
source /home/alberto/.bash_profile

. /home/alberto/programs/miniforge3/bin/activate
conda activate htseq_venv

HTSEQ_SAMPLE_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/pax5_sox8_zfp1_cheng/bowtie2_output/"
HTSEQ_OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/pax5_sox8_zfp1_cheng/bowtie2_htseq_counts/"

echo "HTSeq output folder is $HTSEQ_OUTPUT_FOLDER"
echo "mkdir -p $HTSEQ_OUTPUT_FOLDER"

mkdir -p $HTSEQ_OUTPUT_FOLDER

for f in ${HTSEQ_SAMPLE_FOLDER}/* ; do

    if [ ! -d $f ] ; then
        continue
    fi

    x=${f##*/}
    z=${x//bowtie2_out_}
    echo Starting with sample ${z} ...
    mkdir -p ${HTSEQ_OUTPUT_FOLDER}/htseq_out_${z}

    htseq-count --nonunique all -f bam \
    -r pos --type gene \
    -s no -n 12 \
    ${f}/*.bam \
    $ANNOT \
    > ${HTSEQ_OUTPUT_FOLDER}/htseq_out_${z}/${z}.counts \


    echo "done with sample ${z} ..."

done

echo "Done."

conda deactivate
conda deactivate
