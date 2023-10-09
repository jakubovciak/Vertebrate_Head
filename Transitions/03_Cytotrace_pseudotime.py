#### Markos et al 2023
#### CytoTRACE pseudotime calculation based on integrated data
#### Can be run from bash by `python 03_Cytotrace_pseudotime.py` (with designated conda environment activated!)

# setup ####

import os

print("Current working directory: {0}".format(os.getcwd()))

import scanpy as sc

import numpy as np

import anndata as an

import cellrank as cr

project_name = 'BraFlo100_ABCD'
default_clustering = 'leiden'
annotation='celltype'

# load output of 02_Convert_to_anndata.r
adata=an.read_h5ad('./output/amp_int_mnn_ABCD.h5ad')

# calculate new UMAP for visualisation of pseudotime

sc.pp.neighbors(adata, n_pcs=30, n_neighbors=15,use_rep='corrected')

sc.tl.leiden(adata,resolution=1)

sc.tl.umap(adata,min_dist=0.25,maxiter=500)

# calculate CytoTRACE score and pseudotime
from cellrank.tl.kernels import CytoTRACEKernel

ctk = CytoTRACEKernel(adata,layer='X')

# check umap for consistency of CytoTRACE pseudotime with timepoints assignment
sc.pl.umap(adata,show=False,save='_adata_pt_stage.png',color=['stage','ct_pseudotime'],title=project_name+' UMAP')

sc.pl.violin(adata, keys=["ct_pseudotime"], groupby='stage', rotation=90,show=False,save='_adata_pt.png')

# export metadata to csv table
adata.obs.to_csv('output/adata.obs.csv')




