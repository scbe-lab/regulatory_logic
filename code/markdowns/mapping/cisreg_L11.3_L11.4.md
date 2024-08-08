# Smed cisreg reanalysis: ACME libraries 11.3, 11.4

Libraries 11.3 and 11.4 will be analysed from scratch.

## Concatenating the shallow and deep reads

### L11.3

```sh
cd ~/projects/smed_cisreg/outputs/data/L11.3
# Read1
cat shallow/s11_3_EKDL210003811-1a-10-N501_HCCWVDSX2_L2_1.fq.gz deep/s11_3_EKDL210003811-1a-10-N501_HJNKWDSX2_L3_1.fq.gz > L11_3_1.fq.gz
# Read2
cat shallow/s11_3_EKDL210003811-1a-10-N501_HCCWVDSX2_L2_2.fq.gz deep/s11_3_EKDL210003811-1a-10-N501_HJNKWDSX2_L3_2.fq.gz > L11_3_2.fq.gz
```

### L11.4

```sh
cd ~/projects/smed_cisreg/outputs/data/L11.4
# Read1
cat shallow/s11_4_EKDL210003811-1a-11-N501_HCCWVDSX2_L2_2.fq.gz deep/s11_4_EKDL210003811-1a-11-N501_HJNKWDSX2_L3_2.fq.gz > L11_4_2.fq.gz
# Read2
cat shallow/s11_4_EKDL210003811-1a-11-N501_HCCWVDSX2_L2_2.fq.gz deep/s11_4_EKDL210003811-1a-11-N501_HJNKWDSX2_L3_2.fq.gz > L11_4_2.fq.gz
```

## Cutadapt

### L11.3

```sh
# Read 1
mkdir ~/projects/smed_cisreg/outputs/scrna/L11.3
cd ~/projects/smed_cisreg/outputs/scrna/L11.3
cutadapt -j 4 -m 60 \
-b AGATCGGAAGAG \
-o L11_3_1_trim.fq \
~/projects/smed_cisreg/data/L11.3/L11_3_1.fq.gz \
2> L11_3_1_trim.log 1>&2
# Read 2
cutadapt --report=full -j 4 -m 94 \
--trim-n -b CTGTCTCTTATA ~/projects/smed_cisreg/data/L11.3/L11_3_2.fq.gz | \
cut -c 1-100 | \
grep -B1 -A2 "ATTCG..............$" | gzip -9 > L11_3_2_dephased.fq
```

### L11.4

```sh
# Read 1
mkdir ~/projects/smed_cisreg/outputs/scrna/L11.4
cd ~/projects/smed_cisreg/outputs/scrna/L11.4
cutadapt -j 4 -m 60 \
    -b AGATCGGAAGAG \
    -o L11_4_1_trim.fq \
    ~/projects/smed_cisreg/data/L11.4/L11_4_1.fq.gz \
    2> L11_4_1_trim.log 1>&2
# Read 2
cutadapt --report=full -j 4 -m 94 \
--trim-n -b CTGTCTCTTATA ~/projects/smed_cisreg/data/L11.4/L11_4_2.fq.gz | \
cut -c 1-100 | \
grep -B1 -A2 "ATTCG..............$" | gzip -9 > L11_4_2_dephased.fq
```

## PairFQ

### L11.3

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L11.3
pairfq makepairs \
-f L11_3_1_trim.fq \
-r L11_3_2_dephased.fq.gz \
-fp ./OUT_1_trimmed_pair.fq \
-rp ./OUT_2_trimmed_pair.fq \
-fs ./OUT_1_trimmed_sing.fq \
-rs ./OUT_2_trimmed_sing.fq
```

### L11.4

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L11.4
pairfq makepairs \
-f L11_4_1_trim.fq \
-r L11_4_2_dephased.fq.gz \
-fp ./OUT_1_trimmed_pair.fq \
-rp ./OUT_2_trimmed_pair.fq \
-fs ./OUT_1_trimmed_sing.fq \
-rs ./OUT_2_trimmed_sing.fq
```

## MakeSam

### L11.3

```sh
picard FastqToSam \
F1=OUT_1_trimmed_pair.fq \
F2=OUT_2_trimmed_pair.fq \
O=L11_3_read_pairs.bam \
SM=Smed_L11.3 \
TMP_DIR=./tmp/
```

### L11.4

```sh
picard FastqToSam \
F1=OUT_1_trimmed_pair.fq \
F2=OUT_2_trimmed_pair.fq \
O=L11_4_read_pairs.bam \
SM=Smed_L11.4 \
TMP_DIR=./tmp/
```

## Splitseq pipeline

This will be done using the rink genome so far found at:

`/mnt/sda/alberto/projects/smed_cisreg/data/standard_references`

### L11.3

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L11.3
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/ \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 17000 \
./L11_3_read_pairs.bam
```

### L11.4

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L11.4
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/ \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 18500 \
./L11_4_read_pairs.bam
```
These two go straight up for DGE because here are cells of the HNF knockdown.

## Fishing out only GFP cells

To recover only the GFPi cells, a.k.a. unaffected by the knockdown, we will re-run DigitalExpression using the barcodes of the GFPi cells found at:

```
# L11.3
/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/Sample11/11_4_CellBarcodes_GFP
# L11.4
/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/Sample11/11_3_CellBarcodes_GFP
```

### L11.3

using the barcodes found at:
`/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/Sample11/11_3_CellBarcodes_GFP`

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L11.3_GFPiOnly.txt.gz \
SUMMARY=./L11.3_GFPiOnly_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
CELL_BC_FILE=/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/Sample11/11_3_CellBarcodes_GFP \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```

### L11.4

using the barcodes found at:
`/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/Sample11/11_4_CellBarcodes_GFP`

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L11.4_GFPiOnly.txt.gz \
SUMMARY=./L11.4_GFPiOnly_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
CELL_BC_FILE=/mnt/sda/nkenny/Guo_Assembly/RunSplitseq/Sample11/11_4_CellBarcodes_GFP \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```

## Creating a matrix with all cells (for DGE analysis)


### L11.3

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L11.3_allCells.txt.gz \
SUMMARY=./L11.3_allCells_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp5
```

### L11.4

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L11.4_allCells.txt.gz \
SUMMARY=./L11.4_allCells_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp5
```