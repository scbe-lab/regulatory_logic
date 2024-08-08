import scanpy as sc
import numpy as np
import anndata
import gzip
from scipy.io import mmwrite

# Load the adata file
adata = sc.read_h5ad('/mnt/sda/alberto/projects/smed_cisreg/data/std_ref_cellAnnotation/smed_size_analysis_202306.h5ad')

# Extract the obs and vars
adata.obs.to_csv('/mnt/sda/alberto/projects/smed_cisreg/data/std_ref_cellAnnotation/cells_metadata.csv') # obs
adata.var.to_csv('/mnt/sda/alberto/projects/smed_cisreg/data/std_ref_cellAnnotation/genes_metadata.csv') # var
adata.raw.var.to_csv('/mnt/sda/alberto/projects/smed_cisreg/data/std_ref_cellAnnotation/all_genes_metadata.csv') # var, all

# Get the scaled/normalised counts
scaled_counts = adata.X

# Get the raw counts matrix
raw_counts = adata.raw.X

# Define the output file path
output_file = "/mnt/sda/alberto/projects/smed_cisreg/data/std_ref_cellAnnotation/raw.mtx.gz"

# Save the scaled count matrix as .mtx.gz
with gzip.open(output_file, 'wb') as f:
    mmwrite(f, scaled_counts)

# Save the raw counts matrix as .mtx.gz
with gzip.open(output_file, 'wb') as f:
    mmwrite(f, raw_counts)
