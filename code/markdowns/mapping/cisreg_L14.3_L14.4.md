# Smed cisreg reanalysis: ACME libraries 14.3, 14.4

Libraries 14.3 and 14.4 were pre-paired by Nathan and can be passed directly to the SplitSeq pipeline.

## Splitseq pipeline

This will be done using the rink genome so far found at:

`/mnt/sda/alberto/projects/smed_cisreg/data/standard_references`

### L14.3

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L14.3
ln -s ../../../data/L14.3/read_pairs/Merged7_5_14_3_read_pairs.bam L14_3_read_pairs.bam # this is as of today a .gz file hence why it looks like a broken link in nostromo, probably same for l14_4
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/ \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 6100 \
./L14_3_read_pairs.bam
```
To retain cells with 10 genes minimum, we re-run the DigitalExpression step specifying this particular parameter.

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L14.3.txt.gz \
SUMMARY=./L14.3_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```

### L14.4

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L14.4
ln -s ../../../data/L14.4/read_pairs/Merged7_6_14_4_read_pairs.bam L14_3_read_pairs.bam
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/ \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 6100 \
./L14_4_read_pairs.bam
```
To retain cells with 10 genes minimum, we re-run the DigitalExpression step specifying this particular parameter.

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L14.4.txt.gz \
SUMMARY=./L14.4_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```