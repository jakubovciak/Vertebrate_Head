#### Kozmikova 2023 Amphioxus Head
#### Script used for analysis of 10X matrices
#### Crest markers

# load prerequisities
library(Seurat)
library(patchwork)
library(ggplot2)

# load output of AmphiN5stage.R
AmphiN5stage <- readRDS("../timepoints_rds/Amp_N5.RDS")

# Dorsal Root Ganglia Markers
pl_n5_dot_markers_ganglia<-DotPlot(
  AmphiN5stage,
  features = c("Islet","Trk","SHox2","Tlx","Stac","Slc17A6","Kcna1","Kcnd2","Chrnb2/4","Scrt2","Otof","Hmx","Pou4f3","Prdm12","Skor2","Stmn2","Asic2","Asic1","Plekha8","Lingo2","Neto2","Slc7a14","ScgN","Scg3","Scg5","Npffr2"),
  cols = c('snow2', 'red1'),
  col.min = -2,
  col.max = 10,
  dot.min = 0,
  dot.scale = 10
) + RotatedAxis() + NoLegend() + ggtitle('N5 Dorsal Root Ganglia Markers')

pl_n5_dot_markers_ganglia

# Supplementary Figures

# Neural Crest Markers Expressed Preferentially In Anterior Neural Tube
pl_n5_umap_markers_crest_tube<-FeaturePlot(
  AmphiN5stage,
  features = c("AP2","Snail","FoxD","Dlx","Zic","Otx","Nkx2-2","Fli1","Lhx2/9","Myb","Twist","Six1/2","Insm1","Prdm2","Mef2A","Egr1","Gbx1","Lhx1/5","Mnt","Myc","Tcf7L2","Tcf12"),
  pt.size = 1,
  order = TRUE,
  cols = c('snow2', 'red1'),
  min.cutoff = 0,
  max.cutoff = 2
)

pl_n5_umap_markers_crest_tube

# Neural Crest Markers
pl_n5_dot_markers_crest<-DotPlot(
  AmphiN5stage,
  features = c("Hand","FoxP1","Meis2","SOXD","SOXC","Rxrg","Runx1","RunX1T1","Prox1","Tle","HoxA1","Hmx","Tlx","Kif3A","Nfat","Ebf","Etv1","Slit1","Tnc","Cdh6","Ddx3x","Ddx6","FoxN3","Idh1","Tbx2/3","Maf1","Itgb1","Adam23","Adam17","Adam10","mmp24","FgfA","Akt3","Dusp16","DEusp10","Smad1","Wnt7B","Tcf20","Max","Prdm10","Prdm14","Mitf","Deaf1","Elk1","Ets","FoxO4","IrxA","IrxC","Pbx1","Rarb1","Rarb2","Srf","Tbx1","Tsc22D2","Zeb2","Mycbp2","Vwa5b1"),
  col.min = -2,
  col.max = 10,
  dot.min = 0,
  cols = c('snow2', 'red1'),
  dot.scale = 10
) + RotatedAxis() + NoLegend() + ggtitle('N5 Neural Crest Markers')

pl_n5_dot_markers_crest
