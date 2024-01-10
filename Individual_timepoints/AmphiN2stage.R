#### Markos et al 2023
#### Script used for analysis of 10X matrices
#### Stage N2

# load prerequisities
library(Seurat)
library(patchwork)
library(ggplot2)

# Data load ####
AmphiN2stage.data <-Read10X(data.dir = "../10X_matrices/N2/filtered_feature_bc_matrix/")

# read gene id conversion table
gene_map <-read.csv('../10X_matrices/gene_id_conversion_table.csv', row.names = 1)

# replace BraFlo100 gene id with BraLan3 name where possible
rownames(AmphiN2stage.data) <-
  sapply(rownames(AmphiN2stage.data), USE.NAMES = FALSE, function(loc) {
    if (loc %in% rownames(gene_map)) {
      gene_map[loc, 'Bl3_id_uniq']
    } else {
      loc
    }
  })

# Create seurat object ####
AmphiN2stage <- CreateSeuratObject(
  counts = AmphiN2stage.data,
  project = "AmphiN2stage.ver2",
  min.cells = 1,
  min.features = 1
)

# plot detected genes count vs total count
pl_n2_detected_vs_count <-
  FeatureScatter(AmphiN2stage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('Detected features vs count')
pl_n2_detected_vs_count

# filter cells with too high/low detected genes count
AmphiN2stage <-
  subset(AmphiN2stage, subset = nFeature_RNA > 5 &
           nFeature_RNA < 5000)

# plot again
pl_n2_detected_vs_count_filter <-
  FeatureScatter(AmphiN2stage, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('Detected features vs count, filtered')
pl_n2_detected_vs_count_filter

# Normalize and scale ####
AmphiN2stage <-
  NormalizeData(AmphiN2stage,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

AmphiN2stage <-
  FindVariableFeatures(AmphiN2stage,
                       selection.method = "vst",
                       nfeatures = 5000)

all.genes <- rownames(AmphiN2stage)

AmphiN2stage <- ScaleData(AmphiN2stage, features = all.genes)

# Clustering and dimensionality reduction ####

# run PCA and select PCs for clustering and UMAP
AmphiN2stage <- RunPCA(AmphiN2stage, features = all.genes)
pl_n2_pca<-DimPlot(AmphiN2stage, reduction = "pca") + ggtitle('N2 PCA')
pl_n2_pca

pl_n2_elbow<-ElbowPlot(AmphiN2stage,ndims = 30) + ggtitle('Standard deviation of PCs')
pl_n2_elbow

# run clustering and UMAP
AmphiN2stage <- FindNeighbors(AmphiN2stage, dims = 1:10)
AmphiN2stage <- FindClusters(AmphiN2stage, resolution = 5)
AmphiN2stage <-
  RunUMAP(
    AmphiN2stage,
    dims = 1:10,
    n.neighbors = 10,
    min.dist = 0.6
  )

# inspect clusters and projection
pl_n2_umap1<-DimPlot(
  AmphiN2stage,
  reduction = "umap",
  label = TRUE,
  pt.size = 3,
  order = TRUE
) + NoLegend() + ggtitle('N2 UMAP with clusters')

pl_n2_umap1

# Annotate clusters ####

# assign new cluster names
new.cluster.ids <-
  c(
    "PosteriorSomites",
    "LateralPlateMesoderm",
    "Endoderm",
    "NeuralPlate",
    "Somites",
    "PharyngealEndoderm",
    "Somites",
    "PosteriorSomites",
    "NeuralPlate",
    "AnteriorEndoderm",
    "Endoderm",
    "Notochord",
    "Endoderm",
    "Notochord",
    "Notochord",
    "Notochord",
    "PosteriorSomites",
    "Endoderm",
    "Endoderm",
    "NeuralPlate",
    "Endoderm",
    "Somites",
    "PosteriorSomites",
    "NeuralPlate",
    "PosteriorSomites",
    "AnteriorNeuralPlate",
    "Somites",
    "Unidentified",
    "Somites",
    "PharyngealEndoderm",
    "Endoderm",
    "Somites",
    "NeuralPlate",
    "TailbudMesoderm",
    "PrechordalMesoderm",
    "Endoderm",
    "NeuralPlate",
    "Endoderm",
    "AnteriorNeuralPlate"
  )
names(new.cluster.ids) <- levels(AmphiN2stage)
AmphiN2stage <- RenameIdents(AmphiN2stage, new.cluster.ids)

# inspect new annotation
pl_n2_umap_annot<-DimPlot(
  AmphiN2stage,
  reduction = "umap",
  pt.size = 3,
  order = TRUE
) + ggtitle('N2 UMAP annotated')

pl_n2_umap_annot


# Show selected gene expression for celltypes of interest
AmphiN2stage@active.ident <-
  factor(
    AmphiN2stage@active.ident,
    levels = c(
      "TailbudMesoderm",
      "NeuralPlate",
      "AnteriorNeuralPlate",
      "Notochord",
      "PrechordalMesoderm",
      "PharyngealEndoderm",
      "AnteriorEndoderm",
      "Endoderm",
      "LateralPlateMesoderm",
      "Somites",
      "PosteriorSomites",
      "Unidentified"
    )
  )

pl_n2_dot_gset<-DotPlot(
  AmphiN2stage,
  features = c("Myl12b","Chordin","FoxC","Brachyury1","FoxA1","SoxB1c","Gbx1","Otx","Goosecoid","Nodal","Islet","Bmp2/4","Pax9","Cerberus","Nr6a","Tbx2/3","FoxE4","Rspo2","Axin","Plscr1","Six3/6","Dmbx","Nkx2-1","Nkx2-2","Hex","Six4/5","Six1/2","Fz5/8","Slc6a9","Wnt4","Wnt1","Irxa","Irxb","Fgfa","Fezf","Dkk1/2/4","Dkk3","Lhx2/9","Sim","Shh","Msxlx","Has2","FoxF1","Tbx1","Mnx"),
  cols = c('snow2', 'red1'),
  col.min = -2,
  col.max = 10,
  dot.min = 0,
  dot.scale = 10
) + RotatedAxis()+ NoLegend() + ggtitle('N2 selected genes expression')

pl_n2_dot_gset

# save Figure1 plots
dir.create('../Results/', showWarnings = FALSE)

pdf('../Results/Fig1_N2.pdf',width = 10,height = 7)
pl_n2_umap_annot
pl_n2_dot_gset
dev.off()

# save Seurat object as RDS
dir.create('../timepoints_rds/', showWarnings = FALSE)
saveRDS(AmphiN2stage, file = "../timepoints_rds/Amp_N2.RDS")
