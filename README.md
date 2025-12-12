# regulatory_logic
Code and data for the manuscript "Multimodal single cell analyses reveal gene networks of planarian stem cell differentiation".

## About
This repository hosts the code used to perform the analyses from the manuscript "**Pérez-Posada, A.**; García-Castro, H.; Emili, E.; Guixeras-Fontana, A.; Vanni, V.; Salamanca-Diaz, D.; Arias-Baldrich, C.; van Heeringen, SJ.; Cebrià, F.; Kenny, NJ; Solana, J.. Multimodal single cell analyses reveal gene networks of planarian stem cell differentiation. *Nat Commun*, **16**, 10683 (2025)" , originally titled "The Regulatory Logic of Planarian Stem Cell Differentiation". 

Updated link to Nat Comms: https://www.nature.com/articles/s41467-025-65712-0

DOI: https://doi.org/10.1038/s41467-025-65712-0

Open Access Sharing/Download Link: https://rdcu.be/eR6iZ

Here you can find all the code that was used to generate all the panels from the main and supplementary figures, as well as the majority of supplementary files.

## Data availability:

Main GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274286

- **scRNA-Seq atlas: Seurat objects in .RDS format** (whole atlas and hnf4 knockdown atlas): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274282
- **scATAC-Seq: Seurat object in .RDS format** (together with bigwigs of each cell type): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274281
- bulk RNA-Seq for double knockdown experiments: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299320
- bulk ATAC-Seq: bigwig format https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274280

## Basic structure of the repository
The repository is organised in several folders:

- `code`: the necessary code for the analyses. Within this folder:
  - `code/markdowns`: contains the markdown files recapitulating all the tools and software used from mapping the reads to generation of networks. **Code for the SPLiT-Seq pipeline and ANANSE can be found here**.
  - `code/r_code`: all the R code, including markdowns, static scripts, and functions:
    - `code/r_code/functions`: all the functions used to run the code, **including implementations of other code such as cell normalisation, tau metric, etc; functions to operate with graphs, wrappers for statistical analyses, etc**.
    - `code/r_code/r_scripts`: static scripts such as one-time runs to e.g. extract data from a dataset, re-format files, create some of the supplementary files, etc.
    - `code/r_code/r_markdowns`: **the code for the core analyses**. Although not every markdown needs to be run after the previous one, all markdowns have been ordered numerically following a narrative similar to the one of the manuscript. 
  - `code/scripts`: static scripts used to e.g. deploy the RNA-Seq mapping of the bulk RNA-Seq data (either ours or from the public literature.), running transdecoder, etc.
- `outputs`:
  - `outputs/ananse`: the outputs from ANANSE.
  - `outputs/raw_matrices`: barcodes and matrices, raw.
  - `outputs/bulk_ATAC`: peak calling from the bulk ATAC-Seq
  - `outputs/genomeviz`: tabular data with coordinates for genome visualisation
  - `outputs/gene_annotation`: a table with the eggNOG output and TF annotation; Gene Ontology; and COG categories
  - `outputs/celltype_annotation/`: Tables of cell type annotation
  - `outputs/tables_for_supp/`: Tables that became supplementary files, either with or without further copyediting.
  - `outputs/associate_peaks_genes`: gene/OCR associations
  - `outputs/rda`: R objects that might be of interest for running the markdowns. For example: WGCNA, plotting of ANANSE graphs, TF connectivity analyses...
  - `outputs/hnf_KD`: outputs from the single cell transcriptomics of the HNF knockdown
  - `outputs/motif_analyses`: motif enrichment analyses for promoters of gene modules, OCRs, etc.
  - `outputs/mapped_RNASeq`: mapped (this study) and re-mapped (public) libraries of RNA-Seq data
- `figures`: the .svg and .pdf renders of the figures from the manuscript. Figure legends not included here.
- `data`: static data such as:
  - tables with sample information,
  - tables with data from knockdown in situ hybridisation,
  - data retrieved from querying public databases such as PlanMine
 - `graphics`: legacy directory, mostly used to show some of the featureplots and dimplots that do not appear in the manuscript.


## Notes

Most of these files and the necessary input are available in this repo, but one-click or full-code reproducibility might not be possible: some files have been omitted due to several kind of constraints such as file size limit, intermediate files, or because they are originally from a different publication. The biggest offenders are the .rda and .RDS files that are needed to run the markdowns, and where all of the analyses were generated, as well as the ANANSE networks of the cell types, but also the .bam file from the bulk ATAC-Seq data. **The main objects with the single cell data are available at the GEO associated with this manuscript**, but if other files are needed, feel free to send an email -**they are all available upon request without any restriction**. In that case, please contact Alberto Perez-Posada and/or Jordi Solana, who can provide further information.

Any missing information is fully available upon request without any restriction at ap.posada1[AT]gmail[DOT]com.
