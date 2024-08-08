# Smed cisreg reanalysis: ACME libraries 8.3, 8.4

Starting from pre-trimmed, dephased and paired reads.

## Splitseq pipeline

This will be done using the rink genome so far found at:

`/mnt/sda/alberto/projects/smed_cisreg/data/standard_references`

### L8.3

```sh
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/ \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 19000 \
~/projects/smed_cisreg/data/L8.3/Merged8_3_read_pairs.bam # this is zipped now, see markdown of l14
```

To retain cells with 10 genes minimum, we re-run the DigitalExpression step specifying this particular parameter.

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L8.3.txt.gz \
SUMMARY=./L8.3_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```

### L8.4

```sh
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/ \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 24000 \
~/projects/smed_cisreg/data/L8.4/Merged8_4_read_pairs.bam
```

To retain cells with 10 genes minimum, we re-run the DigitalExpression step specifying this particular parameter.

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L8.4.txt.gz \
SUMMARY=./L8.4_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```