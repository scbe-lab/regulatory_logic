#!/bin/bash
KALLISTO_INDEX="/omics-storage/lab-solana/genomes/Smed/Rink/schMedS3/transdecoder_and_annotation/transdecoder/transcripts.fna.kallisto.index"

which kallisto

# all the single pair-ends
for i in {apc_bcat_brown,bcat_notum_reuter,cdh4_p53_tu,coe_cowles,foxd_vogg,lhx_pitx_currie,mex31_zhu,pax5_sox8_zfp1_unc22_cheng,prmt5_rouhana,yki_lin}; do
    SAMPLE_FOLDER=$(realpath /omics-storage/lab-solana/alberto/DATA/static/reads/${i})
    echo "sample folder is $SAMPLE_FOLDER"

    echo "mkdir $i"
    #cd $i

    KALLISTO_OUTPUT_FOLDER="/omics-storage/lab-solana/alberto/projects/smed_cisreg/outputs/rna/${i}/kallisto_output/raw_kallisto_output"

    echo "output folder is $KALLISTO_OUTPUT_FOLDER"
    echo "mkdir -p $KALLISTO_OUTPUT_FOLDER"
    mkdir -p $KALLISTO_OUTPUT_FOLDER

    FLEN=`grep $i fragment_length.txt | cut -f2`
    FSD=`grep $i fragment_length.txt | cut -f3`
    echo "fragment length for $i is $FLEN"
    echo "fragment sd for $i is $FSD"

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
        -l $FLEN -s $FSD \
        -o ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}
        echo "done with sample ${x} ..."

    done
    echo "Done."

    # cd ../
    pwd

done

