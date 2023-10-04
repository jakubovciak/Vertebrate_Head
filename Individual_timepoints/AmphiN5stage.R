#### Kozmikova 2023 Amphioxus Head
#### Script used for analysis of 10X matrices
#### Stage N5

# load prerequisities
library(Seurat)
library(patchwork)
library(ggplot2)

# Data load ####
AmphiN5stage.data <-Read10X(data.dir = "../10X_matrices/N5/filtered_feature_bc_matrix/")

# read gene id conversion table
gene_map <-read.csv('../10X_matrices/gene_id_conversion_table.csv', row.names = 1)

# replace BraFlo100 gene id with BraLan3 name where possible
rownames(AmphiN5stage.data) <-
  sapply(rownames(AmphiN5stage.data), USE.NAMES = FALSE, function(loc) {
    if (loc %in% rownames(gene_map)) {
      gene_map[loc, 'Bl3_id_uniq']
    } else {
      loc
    }
  })

# Create seurat object ####
AmphiN5stage <- CreateSeuratObject(
  counts = AmphiN5stage.data,
  project = "AmphiN5stage.ver2",
  min.cells = 1,
  min.features = 1
)

# plot detected genes count vs total count
pl_n5_detected_vs_count <-
  FeatureScatter(AmphiN5stage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('Detected features vs count')
pl_n5_detected_vs_count

# filter cells with too high/low detected genes count
AmphiN5stage <-
  subset(AmphiN5stage, subset = nFeature_RNA > 5 &
           nFeature_RNA < 5000)

# plot again
pl_n5_detected_vs_count_filter <-
  FeatureScatter(AmphiN5stage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('Detected features vs count, filtered')
pl_n5_detected_vs_count_filter

# Normalize and scale ####
AmphiN5stage <-
  NormalizeData(AmphiN5stage,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

AmphiN5stage <-
  FindVariableFeatures(AmphiN5stage,
                       selection.method = "vst",
                       nfeatures = 5000)

all.genes <- rownames(AmphiN5stage)

AmphiN5stage <- ScaleData(AmphiN5stage, features = all.genes)

# Clustering and dimensionality reduction ####

# run PCA and select PCs for clustering and UMAP
AmphiN5stage <- RunPCA(AmphiN5stage, features = all.genes)
pl_n5_pca<-DimPlot(AmphiN5stage, reduction = "pca") + ggtitle('N5 PCA')
pl_n5_pca

pl_n5_elbow<-ElbowPlot(AmphiN5stage,ndims = 30) + ggtitle('Standard deviation of PCs')
pl_n5_elbow

# run clustering and UMAP
AmphiN5stage <- FindNeighbors(AmphiN5stage, dims = 1:18)
AmphiN5stage <- FindClusters(AmphiN5stage, resolution = 10)
AmphiN5stage <-
  RunUMAP(
    AmphiN5stage,
    dims = 1:18,
    n.neighbors = 18,
    min.dist = 0.6
  )

# inspect clusters and projection
pl_n5_umap1<-DimPlot(
  AmphiN5stage,
  reduction = "umap",
  label = TRUE,
  pt.size = 3,
  order = TRUE
) + NoLegend() + ggtitle('N5 UMAP with clusters')

pl_n5_umap1

# Annotate clusters ####

# assign new cluster names
new.cluster.ids <-
  c(
    "TailbudMesoderm",
    "LPM1",
    "PharyngealPostDorsalEndo",
    "PharyngealAntEndo1",
    "PharyngealAntEndo2",
    "Somites",
    "Somites",
    "MidGutEndoderm",
    "Unidentified",
    "Notochord",
    "Notochord",
    "PharyngealPostVentralEndo",
    "PharyngealPostVentralEndo",
    "LeftDiverticulum",
    "Somites",
    "NeuralTube",
    "Somites",
    "Notochord",
    "AnteriorNeuralTube",
    "NeuralTube",
    "Somites",
    "MidGutEndoderm",
    "NeuralTube",
    "MidGutEndoderm",
    "HindGutEndoderm",
    "MidGutEndoderm",
    "PharyngealPostVentralEndo",
    "Ectoderm",
    "Somites",
    "HindGutEndoderm",
    "Ectoderm",
    "UnidentifiedPharyngealEndo",
    "NeuralTube",
    "MidGutEndoderm",
    "Somites",
    "MidGutEndoderm",
    "Unidentified",
    "AnteriorNeuralTube",
    "HindGutEndoderm",
    "NeuralCrest-likeCells",
    "Somites",
    "AnteriorTipNotochord",
    "UnidentifiedPharyngealEndo",
    "PharyngealPostDorsalEndo",
    "Notochord",
    "Unidentified",
    "LPM2",
    "PharyngealPostDorsalEndo",
    "HindGutEndoderm",
    "EndostyleEndo",
    "Notochord",
    "Somites"
  )
names(new.cluster.ids) <- levels(AmphiN5stage)
AmphiN5stage <- RenameIdents(AmphiN5stage, new.cluster.ids)

# inspect new annotation
pl_n5_umap_annot <-
  DimPlot(AmphiN5stage, reduction = "umap", pt.size = 3) + ggtitle('N5 UMAP annotated')
pl_n5_umap_annot

# set custom colors
pl_n5_umap_annot_col<-DimPlot(
  AmphiN5stage,
  reduction = "umap",
  pt.size = 3,
  cols = c(
    "TailbudMesoderm" = "orangered1",
    "LPM1" = "deepskyblue",
    "PharyngealPostDorsalEndo" = "lightpink2",
    "PharyngealPostVentralEndo" = "darkgoldenrod",
    "PharyngealAntEndo1" = "orange",
    "PharyngealAntEndo2" = "purple1",
    "Somites" = "green3",
    "MidGutEndoderm" = "palegreen3",
    "Unidentified" = "grey",
    "Notochord" = "royalblue1",
    "LeftDiverticulum" = "turquoise2",
    "NeuralTube" = "mediumorchid1",
    "AnteriorNeuralTube" = "plum1",
    "HindGutEndoderm" = "tan",
    "Ectoderm" = "salmon1",
    "UnidentifiedPharyngealEndo" = "grey40",
    "AnteriorTipNotochord" = "magenta2",
    "NeuralCrest-likeCells" = "skyblue3",
    "LPM2" = "hotpink",
    "EndostyleEndo" = "paleturquoise3",
    "PosteriorNotochord" = "skyblue1"
  )
) + ggtitle('N5 UMAP annotated')

pl_n5_umap_annot_col

# Show selected gene expression for celltypes of interest

AmphiN5stage@active.ident <-
  factor(
    AmphiN5stage@active.ident,
    levels = c(
      "Notochord",
      "AnteriorTipNotochord",
      "AnteriorNeuralTube",
      "NeuralCrest-likeCells",
      "NeuralTube",
      "Somites",
      "TailbudMesoderm",
      "LPM2",
      "LPM1",
      "Ectoderm",
      "Unidentified",
      "HindGutEndoderm",
      "MidGutEndoderm",
      "PharyngealPostDorsalEndo",
      "PharyngealPostVentralEndo",
      "LeftDiverticulum",
      "UnidentifiedPharyngealEndo",
      "EndostyleEndo",
      "PharyngealAntEndo1",
      "PharyngealAntEndo2"
    )
  )

pl_n5_dot_gset<-DotPlot(
  AmphiN5stage,
  features = c(
    "FoxA1",
    "Six1/2",
    "AP2",
    "Hex",
    "Notch",
    "Pax6",
    "Otx",
    "Six4/5",
    "Has2",
    "IrxC",
    "Nkx2-2",
    "Dmbx",
    "IrxB",
    "Sox9",
    "Ptch",
    "Brachyury1",
    "Islet",
    "Gata4/5/6",
    "FoxF1",
    "Blimp",
    "Wnt3",
    "Rspo2",
    "Tbx1",
    "Msxlx",
    "Wnt2",
    "Six3/6",
    "Pitx2",
    "FoxE4",
    "Pax9",
    "FoxQ1",
    "FoxC",
    "Wnt11",
    "Nodal",
    "Shh",
    "Pou3f4",
    "Lhx3",
    "Nkx2-1",
    "Pax2/5",
    "Vent1",
    "Dkk1/2/4",
    "Elav",
    "Fezf"
  ),
  cols = c('snow2', 'red1'),
  col.min = -2,
  col.max = 10,
  dot.min = 0,
  dot.scale = 10
) + RotatedAxis() + NoLegend() + ggtitle('N5 selected genes expression')

pl_n5_dot_gset

# save Seurat object as RDS
dir.create('../timepoints_rds/', showWarnings = FALSE)
saveRDS(AmphiN5stage, file = "../timepoints_rds/Amp_N5.RDS")
