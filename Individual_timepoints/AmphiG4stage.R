#### Kozmikova 2023 Amphioxus Head
#### Script used for analysis of 10X matrices
#### Stage G4

# load prerequisities
library(Seurat)
library(patchwork)
library(ggplot2)

# Data load ####
AmphiG4stage.data <-Read10X(data.dir = "../10X_matrices/G4/filtered_feature_bc_matrix/")

# read gene id conversion table
gene_map <-read.csv('../10X_matrices/gene_id_conversion_table.csv', row.names = 1)

# replace BraFlo100 gene id with BraLan3 name where possible
rownames(AmphiG4stage.data) <-
  sapply(rownames(AmphiG4stage.data), USE.NAMES = FALSE, function(loc) {
    if (loc %in% rownames(gene_map)) {
      gene_map[loc, 'Bl3_id_uniq']
    } else {
      loc
    }
  })

# Create seurat object ####
AmphiG4stage <- CreateSeuratObject(
  counts = AmphiG4stage.data,
  project = "AmphiG4stage.ver2",
  min.cells = 1,
  min.features = 1
)

# plot detected genes count vs total count
pl_g4_detected_vs_count <-
  FeatureScatter(AmphiG4stage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('Detected features vs count')
pl_g4_detected_vs_count

# filter cells with too high/low detected genes count
AmphiG4stage <-
  subset(AmphiG4stage, subset = nFeature_RNA > 5 &
           nFeature_RNA < 5000)

# plot again
pl_g4_detected_vs_count_filter <-
  FeatureScatter(AmphiG4stage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('Detected features vs count, filtered')
pl_g4_detected_vs_count_filter

# Normalize and scale ####
AmphiG4stage <-
  NormalizeData(AmphiG4stage,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

AmphiG4stage <-
  FindVariableFeatures(AmphiG4stage,
                       selection.method = "vst",
                       nfeatures = 5000)

all.genes <- rownames(AmphiG4stage)

AmphiG4stage <- ScaleData(AmphiG4stage, features = all.genes)

# Clustering and dimensionality reduction ####

# run PCA and select PCs for clustering and UMAP
AmphiG4stage <- RunPCA(AmphiG4stage, features = all.genes)

pl_g4_pca<-DimPlot(AmphiG4stage, reduction = "pca") + ggtitle('G4 PCA')
pl_g4_pca

pl_g4_elbow<-ElbowPlot(AmphiG4stage,ndims = 30) + ggtitle('Standard deviation of PCs')
pl_g4_elbow

# run clustering and UMAP
AmphiG4stage <- FindNeighbors(AmphiG4stage, dims = 1:20)
AmphiG4stage <- FindClusters(AmphiG4stage, resolution = 6)

AmphiG4stage <-RunUMAP(
    AmphiG4stage,
    dims = 1:20,
    n.neighbors = 20,
    min.dist = 0.4
  )

# inspect clusters and projection
pl_g4_umap1<-DimPlot(
  AmphiG4stage,
  reduction = "umap",
  label = TRUE,
  pt.size = 3,
  order = TRUE
) + NoLegend() + ggtitle('G4 UMAP with clusters')

pl_g4_umap1

# Annotate clusters ####

# assign new cluster names
new.cluster.ids <-
  c(
    "Ectoderm",
    "NeuralPlate",
    "AnteriorMesendoderm",
    "Ectoderm",
    "NeuralPlate",
    "Ectoderm",
    "AxialMesendoderm",
    "Ectoderm",
    "VentralAnteriorMesendoderm",
    "Ectoderm",
    "LateralMesendoderm",
    "NeuralPlate",
    "Ectoderm",
    "NeuralPlate",
    "Ectoderm",
    "Ectoderm",
    "VentralPosteriorMesendoderm",
    "Ectoderm",
    "Unidentified",
    "PosteriorDorsalMesendoderm",
    "Ectoderm",
    "AxialMesendoderm",
    "Ectoderm",
    "Ectoderm",
    "AxialMesendoderm",
    "Ectoderm",
    "Ectoderm",
    "Ectoderm",
    "PreaxialMesendoderm",
    "LateralMesendoderm",
    "VentralPosteriorMesendoderm",
    "VentralAnteriorMesendoderm",
    "PosteriorDorsalMesendoderm",
    "NeuralPlate",
    "Ectoderm",
    "PosteriorDorsalMesendoderm",
    "NeuralPlate",
    "Ectoderm",
    "PosteriorDorsalMesendoderm",
    "Ectoderm"
  )

names(new.cluster.ids) <- levels(AmphiG4stage)

AmphiG4stage <- RenameIdents(AmphiG4stage, new.cluster.ids)

# inspect new annotation
pl_g4_umap_annot <-
  DimPlot(AmphiG4stage, reduction = "umap", pt.size = 3) + ggtitle('G4 UMAP annotated')
pl_g4_umap_annot

# set custom colors
pl_g4_umap_annot_col<-DimPlot(
  AmphiG4stage,
  reduction = "umap",
  pt.size = 3,
  cols = c(
    "VentralPosteriorMesendoderm" = "orangered1",
    "PosteriorDorsalMesendoderm" = "deepskyblue",
    "AnteriorMesendoderm" = "orange",
    "LateralMesendoderm" = "paleturquoise4",
    "PreaxialMesendoderm" = "maroon1",
    "Unidentified" = "grey",
    "AxialMesendoderm" = "royalblue",
    "NeuralPlate" = "darkorchid1",
    "Ectoderm" = "salmon1",
    "VentralAnteriorMesendoderm" = "palegreen3"
  )
) + ggtitle('G4 UMAP annotated')

pl_g4_umap_annot_col

# Show selected gene expression for celltypes of interest
AmphiG4stage@active.ident <-factor(
    AmphiG4stage@active.ident,
    levels = c(
      "Unidentified",
      "Ectoderm",
      "NeuralPlate",
      "AxialMesendoderm",
      "PreaxialMesendoderm",
      "AnteriorMesendoderm",
      "VentralAnteriorMesendoderm",
      "VentralPosteriorMesendoderm",
      "LateralMesendoderm",
      "PosteriorDorsalMesendoderm"
    )
  )

pl_g4_dot_gset<-DotPlot(
  AmphiG4stage,
  features = c(
    "AP2",
    "Neurogenin",
    "Sox1/2/3",
    "Axin",
    "Fgf8/17/18",
    "Sp5",
    "Delta",
    "Wnt1",
    "Wnt11",
    "Wnt8",
    "Dkk1/2/4",
    "Chordin",
    "Admp",
    "Notch",
    "Bambi",
    "Lefty",
    "Tbx2/3",
    "Lhx3",
    "Blimp",
    "Fz5/8",
    "Tbx15",
    "Hand",
    "Lhx1/5",
    "sFRP2-like",
    "Goosecoid",
    "Vsx2",
    "Six4/5",
    "Otx",
    "FoxA1",
    "Brachyury1",
    "Nodal",
    "Hex",
    "Pitx2",
    "Six3/6",
    "Mrf1",
    "Prdm16",
    "Shh",
    "Vent1"
  ),
  cols = c('snow2', 'red1'),
  col.min = -2,
  col.max = 40,
  dot.min = 0,
  dot.scale = 10
) + RotatedAxis() + ggtitle('G4 selected genes expression')

pl_g4_dot_gset

# save Seurat object as RDS
dir.create('../timepoints_rds/', showWarnings = FALSE)
saveRDS(AmphiG4stage, file = "../timepoints_rds/Amp_G4.RDS")
