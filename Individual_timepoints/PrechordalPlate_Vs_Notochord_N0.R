#### Markos et al 2023
#### Script used for analysis of 10X matrices
#### Prechordal plate vs Notochord markers

# load prerequisities
library(Seurat)
library(patchwork)
library(ggplot2)

# load output of AmphiN0stage.R
AmphiN0stage <- readRDS("../timepoints_rds/Amp_N0.RDS")

DimPlot(AmphiN0stage, reduction = "umap", pt.size = 3)

# Compute differential expression analysis between Prechordal Plate and Notochord and the rest of cells

pp_vs_all_dea <-
  FindMarkers(
    AmphiN0stage,
    ident.1 = "PrechordalPlate",
    ident.2 = NULL,
    only.pos = TRUE,
    logfc.threshold = 0.2,
    min.pct = 0.01,
    max.cells.per.ident = Inf
  )

noto_vs_all_dea <-
  FindMarkers(
    AmphiN0stage,
    ident.1 = "Notochord",
    ident.2 = NULL,
    only.pos = TRUE,
    logfc.threshold = 0.2,
    min.pct = 0.01,
    max.cells.per.ident = Inf
  )

# Subset top significant genes for manual selection of markers
pp_vs_all_dea_sig <- pp_vs_all_dea[pp_vs_all_dea$p_val_adj < 1e-26, ]
noto_vs_all_dea_sig <-
  noto_vs_all_dea[noto_vs_all_dea$p_val_adj < 0.05, ] # We need more relaxed threshold for FDR because of low cell count

### Result tables were manually curated to retain ~65 genes for each celltype meeting designed criteria:
### 1. gene has known function
### 2. gene has down/upregulated expression profile between Notochord and Prechordal Plate

noto_markers<-c("Chmp1a","Gstm3","Creb3l1","Btf3","Eef1b2","Aldh1l2","Minos1","Dpep1","Ube2d2","Adk","Psma3","Atg10","Psmb8","Rnf121","Phb2","Sec62","Ubxn4","Stx6","Tmem65","Imp3","Stub1","Vbp1","Slc35e3","Fam133b","Sssca1","Brap","Faf2","Polr2h","Chp1","Nob1","Trappc5","Bloc1s2","Tomm40","Med19","Sparc","Tomm22","Aimp1","Dohh","Mrpl50","Ndufs4","Ndufb10","Mrpl35","Hccs","Etfa","Mettl10","Timm10","Isca2","Ndufa3","Atp6v1g1","Scube1","Megf6","Ddr2","Scrt1","Wdr83os","Ier3ip1","Kremen","Psat1","Smarcd1","Gamt","Ccdc137","Copg2","Polr2j","Abhd4")

pp_markers<-c("FoxB1","Dpt1","Dpt2","Dpt3","Tinagl1","Angpt2","Pebp1","Nox5","Prkg2","Rhof","Calm1","Actn1","Tenm3","Hmga3","Papln","Pfn1-like","Rspo2","Thbs1","Lamc1","Igfbp7","Tubb4b","Dvl3","Iqgap2","Lrrfip2","Pkdcc","Actnc","Cabyr","Wfikkn1","Gpr84","Dusp7","Togaram1","Actnb","Pcna","Stmn3","Cbx6","Cep152","Hells","Dmd","Ugdh","Cntnap1","Galnt2","Cnr10","Crb1","Tcf21","Cks1b","Rpa3","Crip1","Nuak2","Mcm5","Ankrd28","Mcm6","Rassf8","Ccnd2","Sufu","Reep5","Dnajc17","Bmt2","Coq5","C1qtnf4","Cyp-like","Cyp26a","Cdca7l","C6","Lig1","Kbtbd8l")

# Subset for plotting
subset_N0 <-subset(AmphiN0stage, idents = c("PrechordalPlate", "Notochord"))

# Top Notochord-specific genes

pl_n0_dot_noto <- DoHeatmap(subset_N0,
                            features = noto_markers,
                            size = 4,
                            angle = 0) + NoLegend() + ggtitle('Selected Notochord Markers')

pl_n0_dot_noto

#Top Prechordal Plate specific genes

pl_n0_dot_pp <-DoHeatmap(subset_N0,
                         features = pp_markers,
                         size = 4,
                         angle = 0) + ggtitle('Selected Prechordal Plate Markers')

pl_n0_dot_pp

# Figure S6: Combined heatmap

pl_n0_dot_pp_noto <-DoHeatmap(subset_N0,
                              features = c(noto_markers, pp_markers),
                              size = 4,
                              angle = 0) + ggtitle('Selected Prechordal Plate and Notochord Markers')

pl_n0_dot_pp_noto
