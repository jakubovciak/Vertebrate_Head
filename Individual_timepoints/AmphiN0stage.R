#### Markos et al 2023
#### Script used for analysis of 10X matrices
#### Stage N0

# load prerequisities
library(Seurat)
library(patchwork)
library(ggplot2)

# Data load ####
AmphiN0stage.data <-Read10X(data.dir = "../10X_matrices/N0/filtered_feature_bc_matrix/")

# read gene id conversion table
gene_map <-read.csv('../10X_matrices/gene_id_conversion_table.csv', row.names = 1)

# replace BraFlo100 gene id with BraLan3 name where possible
rownames(AmphiN0stage.data) <-
  sapply(rownames(AmphiN0stage.data), USE.NAMES = FALSE, function(loc) {
    if (loc %in% rownames(gene_map)) {
      gene_map[loc, 'Bl3_id_uniq']
    } else {
      loc
    }
  })

# Create seurat object ####
AmphiN0stage <- CreateSeuratObject(
  counts = AmphiN0stage.data,
  project = "AmphiN0stage.ver2",
  min.cells = 1,
  min.features = 1
)

# plot detected genes count vs total count
pl_n0_detected_vs_count <-
  FeatureScatter(AmphiN0stage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('Detected features vs count')
pl_n0_detected_vs_count

# filter cells with too high/low detected genes count
AmphiN0stage <-
  subset(AmphiN0stage, subset = nFeature_RNA > 5 &
           nFeature_RNA < 5000)

# plot again
pl_n0_detected_vs_count_filter <-
  FeatureScatter(AmphiN0stage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('Detected features vs count, filtered')
pl_n0_detected_vs_count_filter

# Normalize and scale ####
AmphiN0stage <-
  NormalizeData(AmphiN0stage,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

AmphiN0stage <-
  FindVariableFeatures(AmphiN0stage,
                       selection.method = "vst",
                       nfeatures = 5000)

all.genes <- rownames(AmphiN0stage)

AmphiN0stage <- ScaleData(AmphiN0stage, features = all.genes)

# Clustering and dimensionality reduction ####

# run PCA and select PCs for clustering and UMAP
AmphiN0stage <- RunPCA(AmphiN0stage, features = all.genes)

pl_n0_pca<-DimPlot(AmphiN0stage, reduction = "pca") + ggtitle('N0 PCA')
pl_n0_pca

pl_n0_elbow<-ElbowPlot(AmphiN0stage,ndims = 30) + ggtitle('Standard deviation of PCs')
pl_n0_elbow

# run clustering and UMAP
AmphiN0stage <- FindNeighbors(AmphiN0stage, dims = 1:10)
AmphiN0stage <- FindClusters(AmphiN0stage, resolution = 1.5)

AmphiN0stage <-
  RunUMAP(
    AmphiN0stage,
    dims = 1:10,
    n.neighbors = 10,
    min.dist = 0.3
  )

# inspect clusters and projection
pl_n0_umap1<-DimPlot(
  AmphiN0stage,
  reduction = "umap",
  label = TRUE,
  pt.size = 3,
  order = TRUE
) + NoLegend() + ggtitle('N0 UMAP with clusters')

pl_n0_umap1

# Annotate clusters ####

# assign new cluster names
new.cluster.ids <-
  c(
    "AnteriorEndoderm",
    "Ectoderm",
    "Ectoderm",
    "Ectoderm",
    "ParaxialMesoderm",
    "VentralEndoderm",
    "NeuralPlate",
    "Ectoderm",
    "Ectoderm",
    "Ectoderm",
    "ParaxialMesoderm",
    "Ectoderm",
    "LateralEndoderm",
    "PrechordalPlate",
    "Unidentified",
    "TailbudMesoderm",
    "Ectoderm",
    "ParaxialMesoderm",
    "ParaxialMesoderm",
    "NeuralPlate",
    "DorsalBlastoporeLip",
    "Notochord"
  )

names(new.cluster.ids) <- levels(AmphiN0stage)

AmphiN0stage <- RenameIdents(AmphiN0stage, new.cluster.ids)

# inspect new annotation
pl_n0_umap_annot<-DimPlot(AmphiN0stage, reduction = "umap", pt.size = 3) + ggtitle('N0 UMAP annotated')
pl_n0_umap_annot

# set custom colors
pl_n0_umap_annot_col<-DimPlot(
  AmphiN0stage,
  reduction = "umap",
  pt.size = 3,
  cols = c(
    "VentralEndoderm" = "orangered1",
    "ParaxialMesoderm" = "deepskyblue",
    "AnteriorEndoderm" = "orange",
    "LateralEndoderm" = "paleturquoise4",
    "PrechordalPlate" = "maroon1",
    "Unidentified" = "grey",
    "Notochord" = "royalblue",
    "NeuralPlate" = "darkorchid1",
    "Ectoderm" = "salmon1",
    "DorsalBlastoporeLip" = "palegreen3",
    "TailbudMesoderm" = "slategray3"
  )
) + ggtitle('N0 UMAP annotated')

pl_n0_umap_annot_col

# Show selected gene expression for celltypes of interest

AmphiN0stage@active.ident <-
  factor(
    AmphiN0stage@active.ident,
    levels = c(
      "Unidentified",
      "Ectoderm",
      "NeuralPlate",
      "TailbudMesoderm",
      "DorsalBlastoporeLip",
      "Notochord",
      "PrechordalPlate",
      "AnteriorEndoderm",
      "VentralEndoderm",
      "LateralEndoderm",
      "ParaxialMesoderm"
    )
  )


pl_n0_dot_gset<-DotPlot(
  AmphiN0stage,
  features = c(
    "AP2",
    "Snail",
    "Neurogenin",
    "Tcf15",
    "Axin",
    "SoxB1c",
    "Fgf8/17/18",
    "Admp",
    "Nodal",
    "Shh",
    "Six4/5",
    "Chordin",
    "Goosecoid",
    "Dmbx",
    "Six3/6",
    "Six1/2",
    "Fz5/8",
    "Rspo2",
    "Cerberus",
    "Lhx2/9",
    "FoxA1",
    "Brachyury1",
    "Mnx",
    "Dkk1/2/4",
    "Dkk3",
    "Tbx15",
    "Wnt8",
    "Mrf1",
    "Tbx2/3",
    "Hex",
    "Blimp",
    "Prdm16",
    "Nkx2-2",
    "Nkx2-1",
    "Gata4/5/6",
    "IrxC",
    "Wnt1",
    "Wnt4",
    "Wnt5"
  ),
  cols = c('snow2', 'red1'),
  col.min = -2,
  col.max = 10,
  dot.min = 0,
  dot.scale = 10
) + RotatedAxis() + NoLegend() + ggtitle('N0 selected genes expression')

pl_n0_dot_gset

# save Figure1 plots
dir.create('../Results/', showWarnings = FALSE)

pdf('../Results/Fig1_N0.pdf',width = 10,height = 7)
pl_n0_umap_annot_col
pl_n0_dot_gset
dev.off()

# save Seurat object as RDS
dir.create('../timepoints_rds/', showWarnings = FALSE)
saveRDS(AmphiN0stage, file = "../timepoints_rds/Amp_N0.RDS")
