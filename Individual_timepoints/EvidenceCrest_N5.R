#### Markos et al 2023
#### Script used for analysis of 10X matrices
#### Crest markers

# load prerequisities
library(Seurat)
library(patchwork)
library(ggplot2)

# load output of AmphiN5stage.R
AmphiN5stage <- readRDS("../timepoints_rds/Amp_N5.RDS")

AmphiN5stage@active.ident <-
  factor(
    AmphiN5stage@active.ident,
    levels = c(
      "Notochord",
      "AnteriorTipNotochord",
      "AnteriorNeuralTube",
      "NeuralTube",
      "MNCC-likeCells",
      "PMNCC",
      "Somites",
      "TailbudMesoderm",
      "LPM2",
      "LPM1",
      "AntEctoderm",
      "Ectoderm",
      "Unidentified",
      "HindGutEndoderm",
      "MidGut",
      "PharyngealPostDorsalEndo",
      "PharyngealPostVentralEndo",
      "LeftDiverticulum",
      "UnidentifiedPharyngealEndo",
      "Endostyle",
      "PharyngealAntEndo1",
      "PharyngealAntEndo2"
    )
  )

# Figure 3-E
pl_n5_dot_markers_1<-DotPlot(
  AmphiN5stage,
  features = c("Fezf","Otx","FoxD","Lhx1/5","Snail","Dlx","Ap2","Zfhx3","Zic","Draxin","Ebf","Ets","Id2","Idh1","FoxN3","Meis2","Elk1","SoxE","SoxD","SoxC","Zeb2","Btg","Mitf","Slit1","Tle","Rxrg","Wnt7b","Tcf20","Robo1","Cdh6","Arnt2","Atf2","Mef2a","Tmem132c","Ptn","Bmpr1b","Plxna4","Six3/6","Elav","Hand","Tlx","Shox2","Otof","Lingo2","Neto2"),
  cols = c('snow2', 'red1'),
  col.min = -10,
  col.max = 5,
  dot.min = 0,
  dot.scale = 10
) + RotatedAxis() + NoLegend() + ggtitle('Figure 3-E N5 Crest markers')

pl_n5_dot_markers_1

# Figure S7-E
pl_n5_dot_markers_2<-DotPlot(
  AmphiN5stage,
  features = c("Six4/5","Six1/2","Eya","Twist","Pax7a","Pax7b","Msx1","Myc","RhoA/B/C-likea","RhoA/B/C-likeb","Itgb1","Sema6","Max","Mycbp2","Myb","Runx1","Runx1t1","Kif3A","Smad1","Dkk1/2/4","Dkk3","Tcf12","Lrp6-Loc118412841","HoxA1","Adam17","Adam10","Ddx3","Akt3","Mmp24","Rarb2","Deaf1","Maf1","Trk","Chrnb2/4","Prdm12","Hmx","Pdzrn","Gad","Irxc","Kcna1","Prox1"),
  cols = c('snow2', 'red1'),
  col.min = -10,
  col.max = 5,
  dot.min = 0,
  dot.scale = 10
) + RotatedAxis() + NoLegend() + ggtitle('Figure S7-E N5 Crest markers')

pl_n5_dot_markers_2
