# Smed Cisreg: mapping of bulk ATAC-Seq

## Bowtie2

the ATAC_pipe.pl is perl pipeline that uses bowtie2. Works as described in github.com/apposada/ptychodera_cisreg_development/ .

```sh
genome_dir="/mnt/sda/alberto/genomes/Smed/Rink/schMedS3/"

ATAC_dir="/mnt/sda/alberto/colabos/virginia_atac/samples_virginia/"

for i in ${ATAC_dir}/* ; do

  x=${i##*/}
  mkdir -p ${ATAC_dir}/nucf_trim100bp/${x}/

  echo "Treating sample ${x} ..."

  perl ~/programs/ATAC_pipe/ATAC_pipe.pl \
      -f1 ${i}/*_1.fastq.gz -f2 ${i}/*_2.fastq.gz \
      -t 100 \
      -o ${ATAC_dir}/nucf_trim100bp/${x}/${x} \
      -s ${genome_dir}/sizes.genome \
      -i schMedS3_h1.bowtie2 \
      -p 12 -bp $genome_dir \
      -ov_th

  echo "done ${x}"

done

echo "Done."
```

## macs2

```sh
macs2 callpeak -f BED \
               -t ../nucf_trim100bp/atac_vir/atac_vir_nucfree.bed \
               --nomodel \
               --extsize 100 \
               --shift 45 \
               --buffer-size 50000 \
               -g 840213658 \
               -p 0.001 \
               -n atac_vir_macs2peaks \
               --outdir ./macs2_out/

cat macs2_out/atac_vir_macs2peaks_peaks.narrowPeak | \
    mergeBed -d 20 | \
    sortBed | \
    awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"peak_"NR}' \
    > peaks_virginia.bed
```