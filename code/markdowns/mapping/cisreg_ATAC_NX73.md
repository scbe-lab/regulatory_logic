# Smed cisreg reanalysis: 10X scATAC-seq libraries NX73

Starting from scratch, keeping a log of everything that was done.

## Downloading, installing, and preparing the software

```sh
wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1678759737&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Nzg3NTk3Mzd9fX1dfQ__&Signature=Zn2Eil6kQn8kp5PqYDqvKFuL54aYLNLLBgHt4tbcpEvLc8nuRmEAyBNLN8C-5TGSpPH55kkKEdP4aKHrWnsjL4z-9p~ZjGrj8djZ4wsc3CN6hJ8IthFJDY9nY3G9W3a4~09HHAhvnoAzPiL1nLYYvszD382oG6JExDfhlMrI-EuroChNui22AXx0czUX4VE0O--DE3Y7Pdc2mSlCnnHX3ROw40~ASt6xmmBVPkbJ9AKBAJ0yGjlSKK0upyqtD2flarL0uxsJGsgy-JUb36gigFzCgAGh18CSN8qOT1SgExca7lp6l0KJZwuUYhu4A070twREMNzQO9g5vAT4GSZAWg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

md5sum cellranger-atac-2.1.0.tar.gz
tar xzvf cellranger-atac-2.1.0.tar.gz
perl -pi.bak -e "s/PATH=/PATH=~/projects/smed_cisreg/code/software/cellranger/cellranger-atac-2.1.0:/" ~/.bash_profile
source ~/.bash_profile
```

## Preparing the standard references for cellranger-ATAC

### Creating the JASPAR .pfm file

For this, we have associated all Smed genes to huma genes using OrthoFinder. All resolved 1:1 orthologs corresponding to TFs have been properly replaced using 

### Creating consensus peak files

In addition to this, eventually we will map all currently available data as well as novel  bulk ATACseq data to this genome to generate consensus peak files that we can user down the line for generating motif databases as well as high-confidence regions.

[...]

### Creating the config file 
We create a config file that looks like this

```yaml
{
    organism: "Smed_Rink"
    genome: ["Smed_Rink"]
    input_fasta: ["/mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa"]
    input_gtf: ["/mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1_MASKED.gtf"]
}
```

Eventually, we will create a smed database of JASPAR DNA motifs in pfm format, and the yaml will look like this:

```yaml
{
    organism: "Smed_Rink"
    genome: ["Smed_Rink"]
    input_fasta: ["/mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.fa"]
    input_gtf: ["/mnt/sda/alberto/projects/smed_cisreg/data/standard_references/schMedS3_h1.gtf"]
    input_motifs: "/path/to/jaspar/motifs.pfm"
}
```

### Running `cellranger-atac mkref`

After this we create a directory where to put the .cfg and where to run the preparation of the references

```sh
mkdir -p /mnt/sda/alberto/projects/smed_cisreg/data/cellranger_stdrefs
cd /mnt/sda/alberto/projects/smed_cisreg/data/cellranger_stdrefs
#{add here the config file specified above}
cellranger-atac mkref --config=/mnt/sda/alberto/projects/smed_cisreg/data/cellranger_stdrefs/20231107_smed_rink.cfg
```

## Running `cellranger-atac count`

After this we run the cellranger-atac count command that will actually map our reads to the genome.

```sh
cellranger-atac count \
--id=NX73 \
--reference=/mnt/sda/alberto/projects/smed_cisreg/data/cellranger_stdrefs/Smed_Rink/ \
--fastqs=/mnt/sda/alberto/projects/smed_cisreg/data/ATAC/reads/NX73/ \
--localcores=12 \
--localmem=64
```
