# Smed cisreg project: scRNA-seq bam splitting of cell types

## Fishing out reads from cells of interest

```sh
#L01
cat /mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/ACMELib1_merged_mt_MG50_MQ0_DjRem_125above_BC  | awk 'BEGIN {OFS="\t"}{print($1,"smed")}' > l01_smed_barcodes.tsv
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L1/gene_function_tagged.bam) bams/L01_all.bam
samtools index bams/L01_all.bam
sinto filterbarcodes -b bams/L01_all.bam --barcodetag XC -c  l01_smed_barcodes.tsv -p 12 --outdir bams/
mv bams/smed.bam bams/L01_smed.bam

#L02
cat /mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/1_2/ACMELib2_merged_mt_MG50_MQ0_DjRem_125above_BC | awk 'BEGIN {OFS="\t"}{print($1,"smed")}' > l02_smed_barcodes.tsv
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L2/gene_function_tagged.bam) bams/L02_all.bam
samtools index bams/L02_all.bam
sinto filterbarcodes -b bams/L02_all.bam --barcodetag XC -c  l02_smed_barcodes.tsv -p 12 --outdir bams/
mv bams/smed.bam bams/L02_smed.bam

#L03
cat /mnt/sda/nkenny/Guo_Assembly/RunSplitseq/ACME/1_3/ACMELib3_merged_mt_MG50_MQ0_DjRem_125above_BC | awk 'BEGIN {OFS="\t"}{print($1,"smed")}' > l03_smed_barcodes.tsv
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L1_L2_L3/L3/gene_function_tagged.bam) bams/L03_all.bam
samtools index bams/L03_all.bam
sinto filterbarcodes -b bams/L03_all.bam --barcodetag XC -c  l03_smed_barcodes.tsv -p 12 --outdir bams/
mv bams/smed.bam bams/L03_smed.bam

#L08_3
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L8.3/gene_function_tagged.bam) bams/L08_3.bam

#L08_4
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L8.4/gene_function_tagged.bam) bams/L08_4.bam

#L11_3_GFPi
#filter GFPi cells from barcodes
cat /mnt/sda/nkenny/Guo_Assembly/RunSplitseq/Sample11/11_3_CellBarcodes_GFP | awk 'BEGIN {OFS="\t"}{print($1,"gfpi")}' > l11_3_gfpi_barcodes.tsv
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L11.3/gene_function_tagged.bam) bams/L11_3_all.bam
samtools index bams/L11_3_all.bam
sinto filterbarcodes -b bams/L11_3_all.bam --barcodetag XC -c  l11_3_gfpi_barcodes.tsv -p 12 --outdir bams/
mv bams/gfpi.bam bams/L11_3_GFPi.bam

#L11_4_GFPi
#filter GFPi cells from barcodes
cat /mnt/sda/nkenny/Guo_Assembly/RunSplitseq/Sample11/11_4_CellBarcodes_GFP | awk 'BEGIN {OFS="\t"}{print($1,"gfpi")}' > l11_4_gfpi_barcodes.tsv
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L11.4/gene_function_tagged.bam) bams/L11_4_all.bam
samtools index bams/L11_4_all.bam
sinto filterbarcodes -b bams/L11_4_all.bam --barcodetag XC -c  l11_4_gfpi_barcodes.tsv -p 12 --outdir bams/
mv bams/gfpi.bam bams/L11_4_GFPi.bam

#L14_3
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L14.3/gene_function_tagged.bam) bams/L14_3.bam
#L14_4
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L14.4/gene_function_tagged.bam) bams/L14_4.bam
#L23_1
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L23/sublib1/gene_function_tagged.bam) bams/L23_1.bam
#L23_2
ln -s $(realpath ~/projects/smed_cisreg/outputs/scrna/L23/sublib2/gene_function_tagged.bam) bams/L23_2.bam

```

We index these output bam files:

```sh
#all
samtools index bams/L01_smed.bam
samtools index bams/L02_smed.bam 
samtools index bams/L03_smed.bam
samtools index bams/L08_3.bam
samtools index bams/L08_4.bam
samtools index bams/L11_3_GFPi.bam
samtools index bams/L11_4_GFPi.bam
samtools index bams/L14_3.bam
samtools index bams/L14_4.bam
samtools index bams/L23_1.bam
samtools index bams/L23_2.bam
```
## Split all bam files by cells from each cluster (== cell type)

```sh
#separate all by library
mkdir -p outs

for i in {01_smed,02_smed,03_smed,08_3,08_4,11_3_GFPi,11_4_GFPi,14_3,14_4,23_1,23_2} ; do 
    mkdir -p outs/L${i} 
    
    echo "sinto filterbarcodes -b bams/L${i}.bam --barcodetag XC -c smed_cisreg_scrna_cells_cluster_L${i}.tsv --outdir outs/L${i}"

    sinto filterbarcodes -b bams/L${i}.bam \
                         --barcodetag XC \
                         -c smed_cisreg_scrna_cells_cluster_L${i}.tsv \
                         --outdir outs/L${i}
done
```
For each library of the experiment, we have 10-11 bam files corresponding to the reads of each cluster. For each cluster, we can merge the bams of all libraries to retrieve one unique bam file per cell cluster:

```sh
#merge all by celltype
mkdir -p outs/celltypes

for j in {basal_goblet,early_epidermal_progenitors,epidermis,late_epidermal_progenitors,muscle,neoblast,neurons,parenchyma,phagocytes,protonephridia,secretory}; do
    echo "samtools merge -@ 12 - outs/L*/${j}.bam | samtools sort -@ 12 > outs/celltypes/${j}.bam"
    
    samtools merge -@ 12 - outs/L*/${j}.bam | \
    samtools sort -@ 12 > outs/celltypes/${j}.bam
    
    samtools index -@ 12 -b outs/celltypes/${j}.bam
    
    bamCoverage -b outs/celltypes/${j}.bam \
                -o outs/celltypes/${j}.bw
done
```

These bam files are turned into bigwig as well and we can visualise them in R.
