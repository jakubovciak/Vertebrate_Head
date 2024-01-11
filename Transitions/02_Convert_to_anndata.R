#### Markos et al 2023
#### Conversion of integrated SingleCellExperiment object to h5ad for analysis with python
#### Can be run from bash e.g. by `Rscript 02_Convert_to_anndata.R`

# uncomment and set working directory to repository path if run within RStudio:
# setwd('path/to/downloaded/repository')

# Convert integrated data to h5ad format for AnnData input.

library(zellkonverter)
library(Rhdf5lib)
library(HDF5Array)

amp_int <- readRDS('output/amp_merged_sc_int.RDS')

writeH5AD(
  amp_int,
  file = 'output/amp_int_mnn_ABCD.h5ad',
  X_name = 'reconstructed',
  colData = TRUE,
  rowData = FALSE,
  metadata = FALSE
)
