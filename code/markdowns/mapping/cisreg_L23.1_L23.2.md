# Smed cisreg reanalysis: ACME libraries 23

Libraries 23.1 and 23.2 will be analysed from scratch.

## Concatenating the shallow and deep reads

### L23_sublib1

~/projects/smed_cisreg/data/L23/s23_sublib2_1.fq.gz
~/projects/smed_cisreg/data/L23/s23_sublib2_2.fq.gz

```sh
cd ~/projects/smed_cisreg/outputs/data/L23
unzip -j /mnt/sda/NGS_READS/Sample23_size_Sample24_PEG_Jan2022/X204SC21114160-Z01-F001.zip "X204SC21114160-Z01-F001/raw_data/s23_*" -d ./
# Read1
cat s23_1*1.fq.gz > s23_sublib1_1.fq.gz
# Read2
cat s23_1*2.fq.gz > s23_sublib1_2.fq.gz
```

### L23_sublib2

```sh
cd ~/projects/smed_cisreg/outputs/data/L23
# Read1
cat s23_2*1.fq.gz > s23_sublib2_1.fq.gz
# Read2
cat s23_2*2.fq.gz > s23_sublib2_2.fq.gz
```

## Cutadapt

### L23_sublib1

```sh
# Read 1
mkdir ~/projects/smed_cisreg/outputs/scrna/L23/sublib1/
cd ~/projects/smed_cisreg/outputs/scrna/L23/sublib1
cutadapt -j 4 -m 60 \
-b AGATCGGAAGAG \
-o s23_sublib1_1_trim.fq \
~/projects/smed_cisreg/data/L23/s23_sublib1_1.fq.gz \
2> s23_sublib1_1_trim.log 1>&2
# Read 2
cutadapt --report=full -j 4 -m 94 \
--trim-n -b CTGTCTCTTATA ~/projects/smed_cisreg/data/L23/s23_sublib1_2.fq.gz | \
cut -c 1-100 | \
grep -B1 -A2 "ATTCG..............$" | gzip -9 > s23_sublib1_2_dephased.fq
```

### L23_sublib2

```sh
# Read 1
mkdir ~/projects/smed_cisreg/outputs/scrna/L23/sublib2/
cd ~/projects/smed_cisreg/outputs/scrna/L23/sublib2
cutadapt -j 4 -m 60 \
-b AGATCGGAAGAG \
-o s23_sublib1_1_trim.fq \
~/projects/smed_cisreg/data/L23/s23_sublib2_1.fq.gz \
2> s23_sublib2_1_trim.log 1>&2
# Read 2
cutadapt --report=full -j 4 -m 94 \
--trim-n -b CTGTCTCTTATA ~/projects/smed_cisreg/data/L23/s23_sublib2_2.fq.gz | \
cut -c 1-100 | \
grep -B1 -A2 "ATTCG..............$" | gzip -9 > s23_sublib2_2_dephased.fq.gz
```

## PairFQ

### L23_sublib1

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L23/sublib1/
pairfq makepairs \
-f s23_sublib1_1_trim.fq.gz \
-r s23_sublib1_2_dephased.fq.gz \
-fp ./OUT_1_trimmed_pair.fq \
-rp ./OUT_2_trimmed_pair.fq \
-fs ./OUT_1_trimmed_sing.fq \
-rs ./OUT_2_trimmed_sing.fq
```

### L23_sublib2

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L23/sublib2/
pairfq makepairs \
-f s23_sublib2_1_trim.fq.gz \
-r s23_sublib2_2_dephased.fq.gz \
-fp ./OUT_1_trimmed_pair.fq \
-rp ./OUT_2_trimmed_pair.fq \
-fs ./OUT_1_trimmed_sing.fq \
-rs ./OUT_2_trimmed_sing.fq
```

## MakeSam

### L23_sublib1

```sh
picard FastqToSam \
F1=OUT_1_trimmed_pair.fq \
F2=OUT_2_trimmed_pair.fq \
O=L23_sublib1_read_pairs.bam \
SM=Smed_L23sublib1 \
TMP_DIR=./tmp/
```

### L23_sublib2

```sh
picard FastqToSam \
F1=OUT_1_trimmed_pair.fq \
F2=OUT_2_trimmed_pair.fq \
O=L23_sublib2_read_pairs.bam \
SM=Smed_L23sublib2 \
TMP_DIR=./tmp/
```

## Splitseq pipeline

This will be done using the rink genome so far found at:

`/mnt/sda/alberto/projects/smed_cisreg/data/standard_references`

### L23_sublib1

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L23/sublib1/
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/ \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 19000 \
./L23_sublib1_read_pairs.bam
```
To retain cells with 10 genes minimum, we re-run the DigitalExpression step specifying this particular parameter.

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L23.1.txt.gz \
SUMMARY=./L23.1_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```

### L23_sublib2

```sh
cd ~/projects/smed_cisreg/outputs/scrna/L23/sublib2/
Split-seq_pipeline.sh \
-g /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/ \
-r /mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa \
-n 19000 \
./L23_sublib2_read_pairs.bam
```
To retain cells with 10 genes minimum, we re-run the DigitalExpression step specifying this particular parameter.

```sh
dropseq DigitalExpression \
INPUT=./gene_function_tagged.bam \
OUTPUT=./L23.2.txt.gz \
SUMMARY=./L23.2_summary.txt \
LOCUS_FUNCTION_LIST=INTRONIC \
READ_MQ=0 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_DIR=./tmp4
```