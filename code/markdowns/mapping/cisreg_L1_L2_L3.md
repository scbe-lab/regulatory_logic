# Smed cisreg reanalysis: ACME libraries 1,2,3

Starting from the L1,L2,L3. Already trimmed by Nathan:

## Makepairs

### L1

```sh
mkdir ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L1
cd ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L1
pairfq makepairs \
-f ~/projects/smed_cisreg/data/L1_L2_L3/TrimmedReads/S1_1_1_trimmed.fq.gz \
-r  ~/projects/smed_cisreg/data/L1_L2_L3/TrimmedReads/S1_1_2_trimmed.fq.gz \
-fp ./OUT_1_trimmed_pair.fq \
-rp ./OUT_2_trimmed_pair.fq \
-fs ./OUT_1_trimmed_sing.fq \
-rs ./OUT_2_trimmed_sing.fq
```

### L2

```sh
mkdir ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L2/
cd ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L2/
pairfq makepairs \
-f ~/projects/smed_cisreg/data/L1_L2_L3/TrimmedReads/S1_2_1_trimmed_pair.fq.gz \
-r  ~/projects/smed_cisreg/data/L1_L2_L3/TrimmedReads/S1_2_2_trimmed_pair.fq.gz \
-fp ./OUT_1_trimmed_pair.fq \
-rp ./OUT_2_trimmed_pair.fq \
-fs ./OUT_1_trimmed_sing.fq \
-rs ./OUT_2_trimmed_sing.fq
```

### L3

```sh
mkdir ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L3/
cd ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L3/
pairfq makepairs \
-f ~/projects/smed_cisreg/data/L1_L2_L3/TrimmedReads/S1_3_1_trimmed.fq.gz \
-r  ~/projects/smed_cisreg/data/L1_L2_L3/TrimmedReads/S1_3_2_trimmed.fq.gz \
-fp ./OUT_1_trimmed_pair.fq \
-rp ./OUT_2_trimmed_pair.fq \
-fs ./OUT_1_trimmed_sing.fq \
-rs ./OUT_2_trimmed_sing.fq
```

## MakeSam

### L1

```sh
picard FastqToSam \
F1=OUT_1_trimmed_pair.fq \
F2=OUT_2_trimmed_pair.fq \
O=L1_read_pairs.bam \
SM=Smed_L1 \
TMP_DIR=./tmp/
```

### L2

```sh
picard FastqToSam \
F1=OUT_1_trimmed_pair.fq \
F2=OUT_2_trimmed_pair.fq \
O=L2_read_pairs.bam \
SM=Smed_L2 \
TMP_DIR=./tmp/
```

### L3

```sh
picard FastqToSam \
F1=OUT_1_trimmed_pair.fq \
F2=OUT_2_trimmed_pair.fq \
O=L3_read_pairs.bam \
SM=Smed_L2 \
TMP_DIR=./tmp/
```

## Splitseq pipeline

This will be done using the rink genome so far found at:

`/mnt/sda/alberto/projects/smed_cisreg/data/standard_references`

### L1

```sh
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 6300 \
./04_SPLiTseq_Pipeline/L1_read_pairs.bam
```

### L2

```sh
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 6300 \
./04_SPLiTseq_Pipeline/L2_read_pairs.bam
```

### L3

```sh
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 6300 \
./04_SPLiTseq_Pipeline/L3_read_pairs.bam
```

## Digital Expression (keeping only Smed reads)

Although the pipeline.sh returns a DigExpr matrix ready, in this case we must remove the dugesia reads that are also present in these sublibraries. For this, it is better to re-run the digitalExpression step providing a set of whitelisted, *Smichdtea*-only barcodes.

### L1

using the barcodes found at:
`/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/ACMELib1_merged_mt_MG50_MQ0_DjRem_125above_BC`

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L1_SmedOnly.txt.gz \
SUMMARY=./L1_SmedOnly_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
CELL_BC_FILE=/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/ACMELib1_merged_mt_MG50_MQ0_DjRem_125above_BC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```

### L2

Using the barcodes found at:
`/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/1_2/ACMELib2_merged_mt_MG50_MQ0_DjRem_125above_BC`

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L2_SmedOnly.txt.gz \
SUMMARY=./L2_SmedOnly_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
CELL_BC_FILE=/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/1_2/ACMELib2_merged_mt_MG50_MQ0_DjRem_125above_BC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```

### L3

Using the barcodes found at:
`/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/1_3/ACMELib3_merged_mt_MG50_MQ0_DjRem_125above_BC`

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L3_SmedOnly.txt.gz \
SUMMARY=./L3_SmedOnly_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
CELL_BC_FILE=/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/1_3/ACMELib3_merged_mt_MG50_MQ0_DjRem_125above_BC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```