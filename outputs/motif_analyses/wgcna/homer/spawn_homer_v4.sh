#!/bin/bash

GENOME="/mnt/sda/alberto/genomes/Smed/Rink/schMedS3/schMedS3_h1.fa"
PEAKS_DIRECTORY=$(realpath $1)
BG_PEAKS_ALL=$(realpath $2)
OUTDIR=$(realpath $3)
RANDOM_SHUFFLING_SEED=42
NUM_BG_PEAKS=5000

# from https://stackoverflow.com/questions/5914513/shuffling-lines-of-a-file-with-a-fixed-seed
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}


mkdir -p ${OUTDIR}/backgrounds/

for i in ${PEAKS_DIRECTORY}/*.bed ; do
	
	x=${i##*/}
	z=${x%.bed}
	
	#NUM_PEAKS=$(wc -l $i | cut -d " " -f1)
	
	echo "creating bg sequences for $z"
	
	cat $BG_PEAKS_ALL | \
	    grep -v -f $i | \
	    sort --random-source=<(get_seeded_random $RANDOM_SHUFFLING_SEED ) | \
	    head -n $NUM_BG_PEAKS | \
	    sortBed \
	    > ${OUTDIR}/backgrounds/${z}_bg.bed
	
	BG_PEAKS=${OUTDIR}/backgrounds/${z}_bg.bed
	
	echo "treating sample $z"
	
	findMotifsGenome.pl $i $GENOME \
	${OUTDIR}/homer_output_${z} -bg $BG_PEAKS \
	-p 12 -mis 3 -h -mset vertebrates

done
