#### Markos et al 2023
#### Script used for analysis of 10X matrices
#### Zebrafish markers on Amphioxus data

# load prerequisities
library(Seurat)
library(patchwork)
library(ggplot2)

# load output of AmphiG4stage.R

AmphiG4stage<-readRDS("../timepoints_rds/Amp_G4.RDS")

# load output of AmphiN0stage.R
AmphiN0stage<-readRDS("../timepoints_rds/Amp_N0.RDS")

# load output of AmphiN2stage.R
AmphiN2stage<-readRDS("../timepoints_rds/Amp_N2.RDS")

# Subset celltypes of interest
subset_N0 <- subset(AmphiN0stage, idents = c("AnteriorEndoderm","Notochord","PrechordalPlate"))
subset_N2 <- subset(AmphiN2stage, idents = c("PharyngealEndoderm","Notochord","PrechordalMesoderm"))
subset_G4 <- subset(AmphiG4stage,idents = c("AnteriorMesendoderm","AxialMesendoderm","PreaxialMesendoderm"))

subset_N0@active.ident <- factor(subset_N0@active.ident, levels=c("AnteriorEndoderm","PrechordalPlate","Notochord"))
subset_G4@active.ident <- factor(subset_G4@active.ident, levels=c("AnteriorMesendoderm","PreaxialMesendoderm","AxialMesendoderm"))
subset_N2@active.ident <- factor(subset_N2@active.ident, levels=c("PharyngealEndoderm","PrechordalMesoderm","Notochord"))

# Zebrafish Prechordal markers

zebra_prechordal_markers<-c("Chordin","Goosecoid","FoxA1","Id2","Otx","Fz5/8","Dkk1/2/4","Rnd3","sFRP2-like","Six3/6","Blimp","Krt75","Cd63","Nmt2","Bambi","Bmp2/4","Xbp1","Ism1","Timp4","Bhlha15","Ripply2")

pl_g4_dot_zebra_p <-
  DotPlot(
    subset_G4,
    features = zebra_prechordal_markers,
    cols = c('snow2', 'red1'),
    col.min = -2,
    col.max = 10,
    dot.min = 0,
    dot.scale = 5
  ) + NoLegend() + RotatedAxis() + ggtitle('G4 Zebrafish Prechordal Markers')

pl_n0_dot_zebra_p <-
  DotPlot(
    subset_N0,
    features = zebra_prechordal_markers,
    cols = c('snow2', 'red1'),
    col.min = -2,
    col.max = 10,
    dot.min = 0,
    dot.scale = 5
  ) + NoLegend() + RotatedAxis() + ggtitle('N0 Zebrafish Prechordal Markers')

pl_n2_dot_zebra_p <-
  DotPlot(
    subset_N2,
    features = zebra_prechordal_markers,
    cols = c('snow2', 'red1'),
    col.min = -2,
    col.max = 10,
    dot.min = 0,
    dot.scale = 5
  ) + NoLegend() + RotatedAxis() + ggtitle('N2 Zebrafish Prechordal Markers')

pl_g4_dot_zebra_p
pl_n0_dot_zebra_p
pl_n2_dot_zebra_p

# Zebrafish Notochord markers

zebra_noto_markers<-c("Chordin","FoxD","Admp","Id2","FoxA1","Mnx","Brachyury1","Brachyury2","Twist","Cdkn1B/C","Epha","Sall1","Ppp1r14B/C","Shh","Plod1/2","Dag1","P4ha1","Col2A1")

pl_n0_dot_zebra_noto <-
  DotPlot(
    subset_N0,
    features = zebra_noto_markers,
    cols = c('snow2', 'red1'),
    col.min = -2,
    col.max = 10,
    dot.min = 0,
    dot.scale = 5
  ) + NoLegend() + RotatedAxis() + ggtitle('N0 Zebrafish Notochord Markers')

pl_n2_dot_zebra_noto <-
  DotPlot(
    subset_N2,
    features = zebra_noto_markers,
    cols = c('snow2', 'red1'),
    col.min = -2,
    col.max = 10,
    dot.min = 0,
    dot.scale = 5
  ) + NoLegend() + RotatedAxis() + ggtitle('N2 Zebrafish Notochord Markers')

pl_n0_dot_zebra_noto
pl_n2_dot_zebra_noto



