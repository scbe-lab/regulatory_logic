#!/bin/bash
BOWTIE2_INDEX="/omics-storage/lab-solana/genomes/Smed/Rink/schMedS3/schMedS3_h1.bowtie2"
ANNOT="/omics-storage/lab-solana/standard_references/Smed_Rink/Smed_Rink_MASKED.gtf"

# all the single pair-ends
i="pax5_sox8_zfp1_cheng"

SAMPLE_FOLDER=$(realpath /omics-storage/lab-solana/alberto/DATA/static/reads/${i})

echo "sample folder is $SAMPLE_FOLDER"

echo "mkdir $i"

BOWTIE2_OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/${i}/bowtie2_output/"

echo "Bowtie2 output folder is $BOWTIE2_OUTPUT_FOLDER"
echo "mkdir -p $BOWTIE2_OUTPUT_FOLDER"
mkdir -p $BOWTIE2_OUTPUT_FOLDER

FLEN=`grep $i fragment_length.txt | cut -f2`
FSD=`grep $i fragment_length.txt | cut -f3`

echo "fragment length for $i is $FLEN"
echo "fragment sd for $i is $FSD"

echo "BOWTIE2"

for f in ${SAMPLE_FOLDER}/* ; do

    if [ ! -d $f ] ; then
        continue
    fi

    x=${f##*/}
    z=${x}
    echo Starting with sample ${z} ...
    mkdir -p ${BOWTIE2_OUTPUT_FOLDER}/bowtie2_out_${z}

    bowtie2 -x $BOWTIE2_INDEX \
    -r $f/*.gz \
    --quiet -p 12 \
    -X $FLEN | \
    samtools view -b - | samtools sort -m 10G -@ 12 > ${BOWTIE2_OUTPUT_FOLDER}/bowtie2_out_${z}/${z}.bam
    echo "done with sample ${z} ..."

done
    echo "Done."

echo "NOW HTSEQ"

source /home/alberto/.bashrc
source /home/alberto/.bash_profile

. /home/alberto/programs/miniforge3/bin/activate
conda activate htseq_venv

HTSEQ_SAMPLE_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/${i}/bowtie2_output/"
HTSEQ_OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/${i}/bowtie2_htseq_counts/"

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




