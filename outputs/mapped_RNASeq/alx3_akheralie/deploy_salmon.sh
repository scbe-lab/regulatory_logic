#!/bin/bash

# global definitions
INDEX=/omics-storage/lab-solana/genomes/Smed/Rink/schMedS3/transdecoder_and_annotation/transdecoder/transcripts.fna.salmon.index
SAMPLESDIR=/omics-storage/lab-solana/alberto/DATA/static/reads/alx3_akheralie/intact
OUTPUTSDIR=/home/alberto/storage/projects/smed_cisreg/outputs/rna/alx3_akheralie/intact/salmon_output

if [ ! -d $OUTPUTSDIR ] ; then
  echo "global output dir $OUTPUTSDIR does not exist, creating it"
  mkdir -p $OUTPUTSDIR
fi

# loop
for f in ${SAMPLESDIR}/* ; do
  x=${f##*/}
  echo "Starting sample $x ..."

  R1=${f}/*_1.f*q.gz
  R2=${f}/*_2.f*q.gz
  echo "Read 1: $R1"
  echo "Read 2: $R2"

  OUTDIR=${OUTPUTSDIR}/${x}
  echo "OUTDIR is $OUTDIR"
  mkdir -p $OUTDIR

  echo "salmon quant -i $INDEX -l A -1 $R1 -2 $R2 -o $OUTDIR -p 12 "
  salmon quant -i $INDEX -l A -1 $R1 -2 $R2 -o $OUTDIR -p 12
  echo "done sample $x ."

done
echo "Done."
