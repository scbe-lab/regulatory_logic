#!/bin/bash

GENOME="/home/alberto/genomes/Smed/Rink/schMedS3/schMedS3_h1.fa"
PEAKS_DIRECTORY=$(realpath $1)
BG_PEAKS_ALL=$(realpath $2)
OUTDIR=$(realpath $3)

mkdir -p ${OUTDIR}/backgrounds/

for i in ${PEAKS_DIRECTORY}/*.bed ; do
	
	x=${i##*/}
	z=${x%.bed}
	
	#echo "creating bg sequences for $z"
	
	#cat $BG_PEAKS_ALL | \
	#    sortBed | \
	#    grep -v -f $i \
	#    > ${OUTDIR}/backgrounds/${z}_bg.bed
	
	#BG_PEAKS=${OUTDIR}/backgrounds/${z}_bg.bed
	
	echo "treating sample $z"
	
	findMotifsGenome.pl $i $GENOME \
	${OUTDIR}/homer_output_${z} -bg $BG_PEAKS_ALL \
	-p 12 -mset vertebrates -size 200

done

