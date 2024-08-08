# Smed Cisreg: ANANSE

## Preparation

### Installing genomes

```sh
activate_miniconda3_venv
mamba create -n ananse_venv -c bioconda -c conda-forge ananse #agat too?
conda activate ananse_venv
genomepy install -a GRCh38.p13 -p Ensembl
genomepy install -a GRCm38.p6 -p Ensembl
genomepy install -a BraLan2 -p Ensembl # or ncbi?
# or, copy these genome installations from a previous installation
```

### Preparing the Smed genome

```sh
cd ~/.local/share/genomes/
mkdir -p Smed
cd ./Smed/
# maybe none of this is needed:
# cp ~/../Standard_References/Smed_Rink/* .
# mv schMedS3_h1.fa Smed.fa
# samtools faidx -i Smed.fa -o Smed.fa.fai
# cut -f1,2 Smed.fa.fai > Smed.fa.sizes
# mv schMedS3_h1.gtf Smed.annotation.gtf
# agat_convert_sp_gff2bed.pl --gff Smed.annotation.gtf -o Smed.annotation.bed
cp ${PATH_TO_GENOME} ./Smed.fa # line 112 in genomepy/genome/__init__.py _parse_filename requires a .fa to be present in the directory. hopefully just the .fa?
awk 'BEGIN {OFS="\t"} {print(">"$2,">"$1)}' ${PATH_TO_ROSETTA} > transcript_gene.dct
python2 ~/programs/dictionary_generator/replace_name_fasta.py -i ${PATH_TO_PROTEOME} -l transcript_gene.dct -o pep_seqs_translated #translate from transcript to protein ID
cat pep_seqs_translated | perl -pe "s/\>([A-Za-z0-9]+)/\>\1|\1/" > Smed.pep.fa
rm pep_seqs_translated transcript_gene.dct
# get the different proteomes in the .local/share/genomes directory
```

### Running Motif2Factors

see gimmemotifs.readthedocs.io

We used the following species:

```
(((((((((GRCh38.p13,GRCm38.p6),(DANRE,Locu)),BraLan2),STRPU),((((((Smed,Spol),Djap),Bsem),Avag),(Ofus,Cgig)),(Rvar,DROME))),Hmia),NEMVE),Aque),Cowc);
```

And this was the command. Since we also have JASPAR matrices from having used the inference tool (ref. JASPAR 2024, ref. Similarity regression), we will be using the JASPAR database instead of the gimemmotifs database.

```sh
gimme motif2factors --new-reference Smed --database JASPAR2020 --lenient --tmpdir ./tmp_m2f/ --ortholog-reference Avag Cgig Rvar Hmia Aque Cowc Djap BraLan2 DANRE Locu STRPU Spol Bsem Ofus DROME NEMVE --threads 12 --outdir ./smed_motif2factors_more_refs_JASPAR_2020
```

### Running inference tool from JASPAR

```sh
cd ~/projects/smed_cisreg/outputs/jaspar_inference
mkdir ./jaspar_inference
cd ./jaspar_inference
fastaliner ~/.local/share/genomes/Smed/Smed.pep.fa > all.pep 
grep -A1 -f tfs_all_names.txt all.pep | perl -pe "s/\-\-\n//" > tfs.pep
git clone https://github.com/wassermanlab/JASPAR-inference-tool.git
activate_miniconda3_venv
mamba env create -n inference_venv -f JASPAR-inference-tool/conda/environment.yml 
conda activate inference_venv
python ./JASPAR-inference-tool/infer_profile.py --latest tfs.pep > inferred_matrix.tsv
```

### Merging motif2factors, profile inference from JASPAR, and Neiro et al. 2022 Motif evidence

(see R)

```r
# Load DB of all jaspar motifs used by ANANSE
ananse_jaspar2020_db <- read.delim2(
  "/mnt/sda/alberto/programs/miniconda3/envs/ananse_venv/lib/python3.10/site-packages/data/motif_databases/JASPAR2020.motif2factors.txt",
  sep = "\t",
  header = TRUE
)
# ANANSE has the motif and the TF ids joined in a single column separated by "_"
# Add a column of just the motif ID
ananse_jaspar2020_db$MAid <- gsub("_.*","",ananse_jaspar2020_db$Motif)

# Load motif2factors output
smed_m2f <-read.delim2("~/projects/smed_cisreg/outputs/ananse/m2f/smed_motif2factors_more_refs_JASPAR_2020/Smed.JASPAR2020.motif2factors.txt")

# some gene IDs are carrying over the ">" symbol, for some reason
smed_m2f$Factor <- sub(">","",smed_m2f$Factor)

# Load profile inference results
smed_proInf_results <- read.delim2("~/projects/smed_cisreg/outputs/jaspar_inference/inferred_matrix.tsv")

# we remove dupped ids separated by "|"
smed_proInf_results$Query <- gsub("\\|.*","",smed_proInf_results$Query)

# we merge with the jaspar db so we can give these motif predictions a matrix_TF id which will facilitate the joining of the tables
smed_proInf_results_m2f <-
  merge(
    smed_proInf_results,
    ananse_jaspar2020_db,
    by.x = "TF.Matrix",
    by.y = "MAid",
    all.x = TRUE
  )

# profile inference was done using JASPAR v2024 matrices and m2f is from 2020. Some matrices are not in the m2f file. We miss 5 TFs this way. But it should be ok because we are also doing gimme's m2f.
# for these tfs, "Y" means "taken from profile inference tool"
smed_inference <- smed_proInf_results_m2f[complete.cases(smed_proInf_results_m2f),][,c(8,2,10,11)] 

# Give it the same column names as m2f output to allow joining the tables
colnames(smed_inference) <- colnames(smed_m2f)

smed_m2f_merged <-
  rbind(
    smed_m2f,
    smed_inference
  )

# which of these motifs have now found an "ortholog"? (predicted by JASPAR's inference). We'll have to remove these entries
filt_out <- 
  which(
    smed_m2f_merged$Motif %in% smed_m2f_merged$Motif[smed_m2f_merged$Factor != "NO ORTHOLOGS FOUND"] & # which row has a motif with an assigned ortholog
      smed_m2f_merged$Factor == "NO ORTHOLOGS FOUND" # and yet it says in this row "NO ORTHOLOGS FOUND"
  )

smed_m2f_merged <- smed_m2f_merged[-filt_out,]

# Load the tfs with an assigned JASPAR matrix by Neiro et al., 2022
smed_tfsneiro <- tfs_all[!(tfs_all$Jaspar.MatrixID %in% c(".","-", "0")),c(1,21)]

# we merge with the jaspar db so we can, again, give them a matrix_TF id which will facilitate the joining of the tables
smed_tfsneiro_m2f <- 
  merge(
    smed_tfsneiro,
    ananse_jaspar2020_db,
    by.x = 2,
    by.y = 5,
    all.x = TRUE
  )

smed_tfsneiro_m2f <- smed_tfsneiro_m2f[complete.cases(smed_tfsneiro_m2f),] # 201 out of 286... these motif matrices also come from JASPAR 2020 supposedly. Where are the missing matrices?

# again we rearrange and tidy to do the joining of tables
smed_tfsneiro_m2f <- smed_tfsneiro_m2f[,c(3,2,5,6)]
colnames(smed_tfsneiro_m2f) <- colnames(ananse_jaspar2020_db)[1:4]

smed_m2f_merged <-
  rbind(
    smed_m2f_merged,
    smed_tfsneiro
  )
# which of these motifs have now found an "ortholog"? (predicted by Neiro et al. when they used JASPAR profile inference)
filt_out <- 
  which(
    smed_m2f_merged$Motif %in% smed_m2f_merged$Motif[smed_m2f_merged$Factor != "NO ORTHOLOGS FOUND"] & # which row has a motif with an assigned ortholog
      smed_m2f_merged$Factor == "NO ORTHOLOGS FOUND" # and yet it says in this row "NO ORTHOLOGS FOUND"
  )

smed_m2f_merged <- smed_m2f_merged[-filt_out,]

# load motif2factors from a previous run to do some comparisons
m2f_previous <- read.delim2("~/projects/smed_cisreg/outputs/ananse/m2f/_obs/smed_motif2factors/Smed.gimme.vertebrate.v5.0.motif2factors.txt")
length(which(unique(m2f_previous$Factor) %in% tfs_all$gene)) / nrow(tfs_all)

length(which(unique(smed_m2f_merged$Factor) %in% tfs_all$gene)) # 401 out of 665. Previously we had 196 out of 665
length(which(unique(smed_m2f_merged$Factor) %in% tfs_all$gene)) / nrow(tfs_all)

# should I just keep the ones we got in the curated TF list or should I just add all? I think I should keep the ones we have in our curated list. Should I have used the same roster of species in the TFanimalDB evidence as here with ANANSE?

write.table(
  smed_m2f_merged,
  "outputs/ananse/m2f/smed_motif2factors_more_refs_JASPAR_2020/merge_m2f_inference_neiroetal/Smed.JASPAR2020.motif2factors.txt",
  sep = "\t",
  quote = F,
  row.names = F
)
```

### get peak summits & max read depths

These cells below are an extct copy of the .ipynb located at `code/markdowns/atac_counts.ipynb`, which in itself is an updated running copy of the one from Siebren back in mid 2022.

```python
import pandas as pd
import numpy as np
import subprocess as sp
import os	
```

The peak file in the raw data has two issues:
 1. the peaks have variable widths
 2. some "peaks" have very few reads

Here we make a pileup of all BAM reads to find the number of reads, and the summit of reads under each peak.

Next, we filter the peakset to peaks which contain at least four reads in total.

Finally, we normalize the peak width too 200bp, centered on the summit.

```python
peakfile = "/mnt/sda/alberto/projects/smed_cisreg/outputs/scatac/NX73/outs/raw_peak_bc_matrix/peaks.bed"
df = pd.read_table(peakfile, comment="#", header=None)
df.columns = ["chrom", "start", "end"]
df
```

```python
# no comments/headers
df.to_csv("/mnt/sda/alberto/projects/smed_cisreg/outputs/ananse/prep/peaks.bed", sep="\t", index=False, header=False)
```

```sh
bedtools slop \
-b 1000 \
-i ~/projects/smed_cisreg/outputs/ananse/prep/peaks.bed \
-g ~/genomes/Smed/Rink/schMedS3/sizes.genome \
> ~/projects/smed_cisreg/outputs/ananse/prep/peaks_extended.bed
```

```sh
# -B: Bam already filtered
# -A: Bam already filtered
samtools mpileup \
-A -B \
-l ~/projects/smed_cisreg/outputs/ananse/prep/peaks_extended.bed \
-o ~/projects/smed_cisreg/outputs/ananse/prep/peak_extended_pileup.tsv \
~/projects/smed_cisreg/outputs/ananse/prep/bams/possorted_bam.bam
```

```python
peakfile = "/mnt/sda/alberto/projects/smed_cisreg/outputs/scatac/NX73/outs/raw_peak_bc_matrix/peaks_extended_noname.bed"
df = pd.read_table(peakfile, comment="#", header=None)
pu = pd.read_table(
    "~/projects/smed_cisreg/outputs/ananse/prep/peak_extended_pileup.tsv", 
    header=None,
    usecols=[0,1,3],
    names=["chrom", "pos", "count"],
)
```

```python
summits = []
depths = []
last_chrom = ""
for idx, (chrom, start, end) in df.iterrows():
    if last_chrom != chrom:
        print("chrom:", chrom)  # print progress
        chrom_pu = pu[pu["chrom"] == chrom]
        last_chrom = chrom
    
    sub_pu = chrom_pu[chrom_pu["pos"].between(start, end)]
    maxima = sub_pu["count"] == sub_pu["count"].max()
    if sum(maxima) == 0:
        summit = None
    elif sum(maxima) == 1:
        summit = int(sub_pu[maxima]["pos"])
    else:
        idx = sub_pu[maxima].index
        middle_summit = idx[len(idx) // 2]
        summit = int(sub_pu.loc[middle_summit]["pos"])
    depth = sub_pu["count"].max()
    
    depths.append(depth)
    summits.append(summit)
```

```python
df2 = df.copy()
df2["summit"] = summits
df2["depth"] = depths
df2 = df2.dropna()
df2["summit"] = df2["summit"].astype(np.int64)
df2["depth"] = df2["depth"].astype(np.int64)
```

```python
df2.to_csv("~/projects/smed_cisreg/outputs/ananse/prep/peak_extended_summits.tsv", index=False, sep="\t")
```

Visualise peaks per read depth
```python
df2 = pd.read_csv("~/projects/smed_cisreg/outputs/ananse/prep/peak_extended_summits.tsv", sep="\t")

depth = list(range(0, df2["depth"].max()))
n_regions = []
for n in depth:
    nr = sum(df2['depth'] == n)
    n_regions.append(nr)
    if n <= 10:
        print(
            f"regions with depth {n}: {nr}"
        )
```

```python
import matplotlib
%matplotlib inline
```

```python
_ = pd.DataFrame({"reads": depth, "number of reads under peaks": n_regions}).plot(x="reads", xlim=(0, 50))#, ylim=(0, 20_000))
```

```python
min_depth = list(range(0, df2["depth"].max()))
n_regions = []
for n in min_depth:
    nr = sum(df2['depth'] >= n)
    n_regions.append(nr)
```

```python
_ = pd.DataFrame({"min_depth": min_depth, "remaining_regions": n_regions}).plot(x="min_depth", xlim=(0, 50))#, ylim=(0, 20_000))
```

Filter peaks by min read depth
```python
min_reads = 5 # should be 4??
```

```python
df3 = df2.copy()
df3 = df3[df3['depth'] >= min_reads]
# standardize the peak widths
df3["start"] = df3["summit"]-100
df3["end"] = df3["summit"]+100
df3 = df3[["chrom", "start", "end"]]
df3
```

```python
df3.to_csv(f"~/projects/smed_cisreg/outputs/ananse/prep/peaks_normalized_min{min_reads}.bed", index=False, header=False, sep="\t")
```

## Running ANANSE standalone (one by one)

### Preparation of RNA and ATAC counts for running separate runs of ANANSE per cell type

We will use the bam and the summit-centered peak file. We will let ANANSE parse the bam.

We do the same with the RNA counts from the pseudobulk table. We need to update this!!:

```r
# in R during markdown 08b:
smed_cisreg_scrna_pseudobulk_broad <- 
  pseudobulk(
    x = smed_cisreg_scrna@assays$RNA@counts,
    identities = broad_idents
    )

for (cell in colnames(smed_cisreg_scrna_pseudobulk_broad)){
  write.table(
    smed_cisreg_scrna_pseudobulk_broad[,cell],
    file = paste0("outputs/ananse/data/rna_",cell,".tsv"),
    row.names = TRUE, col.names = FALSE, sep = "\t"
    )
  message("Done ",cell)
  }
```

### ANANSE binding

```sh
cd ~/projects/smed_cisreg/outputs/ananse/
activate_miniconda3_venv
conda activate ananse_venv

M2F_DB="m2f/smed_motif2factors_more_refs_JASPAR_2020/merge_m2f_inference_neiroetal/core_common_tfs/Smed.JASPAR2020.pfm"
GENOME="~/.local/share/genomes/Smed/Smed.fa"
ALL_PEAKS="~/projects/smed_cisreg/outputs/ananse/prep/peaks_normalized_min5.bed"
# originally I used : ~/projects/smed_cisreg/outputs/scatac/NX73/outs/raw_peak_bc_matrix/peaks.bed
# maybe use this one: ~/projects/smed_cisreg/outputs/ananse/prep/peak_extended_summits.bed
# maybe use this one: ~/projects/smed_cisreg/outputs/ananse/prep/peaks_normalized_min5.bed

# neoblast
ananse binding -A data/atac/bams/neoblast.bam -g $GENOME -p $M2F_DB -o outs/binding/neoblast -r $ALL_PEAKS -n 12
# eep
ananse binding -A data/atac/bams/early_epidermal_progenitors.bam -g $GENOME -p $M2F_DB -o outs/binding/eep -r $ALL_PEAKS -n 12
# lep
ananse binding -A data/atac/bams/late_epidermal_progenitors.bam -g $GENOME -p $M2F_DB -o outs/binding/lep -r $ALL_PEAKS -n 12
# epidermis
ananse binding -A data/atac/bams/epidermis.bam -g $GENOME -p $M2F_DB -o outs/binding/epidermis -r $ALL_PEAKS -n 12
# phagocytes
ananse binding -A data/atac/bams/phagocytes.bam -g $GENOME -p $M2F_DB -o outs/binding/phagocytes -r $ALL_PEAKS -n 12
# basalgoblet
ananse binding -A data/atac/bams/basal_goblet.bam -g $GENOME -p $M2F_DB -o outs/binding/basalgoblet -r $ALL_PEAKS -n 12
# muscle
ananse binding -A data/atac/bams/muscle.bam -g $GENOME -p $M2F_DB -o outs/binding/muscle -r $ALL_PEAKS -n 12
# neuron
ananse binding -A data/atac/bams/neurons.bam -g $GENOME -p $M2F_DB -o outs/binding/neuron -r $ALL_PEAKS -n 12
# parenchyma
ananse binding -A data/atac/bams/parenchyma.bam -g $GENOME -p $M2F_DB -o outs/binding/parenchyma -r $ALL_PEAKS -n 12
# protonephridia
ananse binding -A data/atac/bams/protonephridia.bam -g $GENOME -p $M2F_DB -o outs/binding/protonephridia -r $ALL_PEAKS -n 12
# secretory
ananse binding -A data/atac/bams/secretory.bam -g $GENOME -p $M2F_DB -o outs/binding/secretory -r $ALL_PEAKS -n 12
```

### ANANSE network

```sh
GENOME="~/.local/share/genomes/Smed/Smed.fa"
ANNOTATION="~/.local/share/genomes/Smed/Smed.annotation.gtf"
ALL_PEAKS="~/projects/smed_cisreg/outputs/ananse/prep/peaks_normalized_min5.bed"
# neoblast
ananse network -n 24 -e data/rna_neoblast.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/neoblast.network outs/binding/neoblast/binding.h5
# eep
ananse network -n 24 -e data/rna_eep.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/eep.network outs/binding/eep/binding.h5
# lep
ananse network -n 24 -e data/rna_lep.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/lep.network outs/binding/lep/binding.h5
# epidermis
ananse network -n 24 -e data/rna_epidermis.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/epidermis.network outs/binding/epidermis/binding.h5
# phagocytes
ananse network -n 24 -e data/rna_phagocytes.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/phagocytes.network outs/binding/phagocytes/binding.h5
# basalgoblet
ananse network -n 24 -e data/rna_BG.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/basalgoblet.network outs/binding/basalgoblet/binding.h5
# muscle
ananse network -n 24 -e data/rna_muscle.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/muscle.network outs/binding/muscle/binding.h5
# neuron
ananse network -n 24 -e data/rna_neuron.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/neuron.network outs/binding/neuron/binding.h5
# parenchyma
ananse network -n 24 -e data/rna_parenchyma.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/parenchyma.network outs/binding/parenchyma/binding.h5
# protonephridia
ananse network -n 24 -e data/rna_protonephridia.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/protonephridia.network outs/binding/protonephridia/binding.h5
# secretory
ananse network -n 24 -e data/rna_secretory.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/secretory.network outs/binding/secretory/binding.h5
```

### ANANSE influence

(here goes the differential gene expression analysis between cell types; see markdown 08c)
```r
# (here goes the differential gene expression analysis between cell types; see markdown 08c)
```

```sh
# eep_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/eep.network -d outs/influence/diff_eep_neoblast.tsv -o outs/influence/infl_eep_neoblast.txt -n 1 -i 100000
# lep_eep
#ananse influence -s outs/network/eep.network -t outs/network/lep.network -d outs/influence/diff_lep_eep.tsv -o outs/influence/infl_lep_eep.txt -n 1 -i 100000
# epidermis_lep
#ananse influence -s outs/network/lep.network -t outs/network/epidermis.network -d outs/influence/diff_epidermis_lep.tsv -o outs/influence/infl_epidermis_lep.txt -n 1 -i 100000
# epidermis_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/epidermis.network -d outs/influence/diff_epidermis_neoblast.tsv -o outs/influence/infl_epidermis_neoblast.txt -n 1 -i 100000
# phagocytes_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/phagocytes.network -d outs/influence/diff_phagocytes_neoblast.tsv -o outs/influence/infl_phagocytes_neoblast.txt -n 1 -i 100000
# basalgoblet_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/basalgoblet.network -d outs/influence/diff_BG_neoblast.tsv -o outs/influence/infl_basalgoblet_neoblast.txt -n 1 -i 100000
# muscle_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/muscle.network -d outs/influence/diff_muscle_neoblast.tsv -o outs/influence/infl_muscle_neoblast.txt -n 1 -i 100000
# neuron_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/neuron.network -d outs/influence/diff_neuron_neoblast.tsv -o outs/influence/infl_neuron_neoblast.txt -n 1 -i 100000
# parenchyma_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/parenchyma.network -d outs/influence/diff_parenchyma_neoblast.tsv -o outs/influence/infl_parenchyma_neoblast.txt -n 1 -i 100000
# protonephridia_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/protonephridia.network -d outs/influence/diff_protonephridia_neoblast.tsv -o outs/influence/infl_protonephridia_neoblast.txt -n 1 -i 100000
# secretory_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/secretory.network -d outs/influence/diff_secretory_neoblast.tsv -o outs/influence/infl_secretory_neoblast.txt -n 1 -i 100000
```


```sh
mkdir -p outs/influence/250k/
# eep_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/eep.network -d outs/influence/diff_eep_neoblast.tsv -o outs/influence/250k/infl_eep_neoblast_250k.txt -n 1 -i 250000
# lep_eep
#ananse influence -s outs/network/eep.network -t outs/network/lep.network -d outs/influence/diff_lep_eep.tsv -o outs/influence/250k/infl_lep_eep_250k.txt -n 12 -i 250000
# epidermis_lep
#ananse influence -s outs/network/lep.network -t outs/network/epidermis.network -d outs/influence/diff_epidermis_lep.tsv -o outs/influence/250k/infl_epidermis_lep_250k.txt -n 12 -i 250000
# epidermis_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/epidermis.network -d outs/influence/diff_epidermis_neoblast.tsv -o outs/influence/250k/infl_epidermis_neoblast_250k.txt -n 12 -i 250000
# phagocytes_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/phagocytes.network -d outs/influence/diff_phagocytes_neoblast.tsv -o outs/influence/250k/infl_phagocytes_neoblast_250k.txt -n 12 -i 250000
# basalgoblet_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/basalgoblet.network -d outs/influence/diff_BG_neoblast.tsv -o outs/influence/250k/infl_basalgoblet_neoblast_250k.txt -n 12 -i 250000
# muscle_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/muscle.network -d outs/influence/diff_muscle_neoblast.tsv -o outs/influence/250k/infl_muscle_neoblast_250k.txt -n 12 -i 250000
# neuron_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/neuron.network -d outs/influence/diff_neuron_neoblast.tsv -o outs/influence/250k/infl_neuron_neoblast_250k.txt -n 12 -i 250000
# parenchyma_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/parenchyma.network -d outs/influence/diff_parenchyma_neoblast.tsv -o outs/influence/250k/infl_parenchyma_neoblast_250k.txt -n 12 -i 250000
# protonephridia_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/protonephridia.network -d outs/influence/diff_protonephridia_neoblast.tsv -o outs/influence/250k/infl_protonephridia_neoblast_250k.txt -n 12 -i 250000
# secretory_neoblast
ananse influence -s outs/network/neoblast.network -t outs/network/secretory.network -d outs/influence/diff_secretory_neoblast.tsv -o outs/influence/250k/infl_secretory_neoblast_250k.txt -n 12 -i 250000
```
