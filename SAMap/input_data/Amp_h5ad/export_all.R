library(Seurat)
library(zellkonverter)
library(batchelor)
library(scater)
library(dplyr)
library(Rhdf5lib)
library(HDF5Array)

coarse_annot<-read.csv('amp_coarse_annot.csv')

amp_list <- list(
  G4 = readRDS('../../../Vertebrate_Head/timepoints_rds/Amp_G4.RDS'),
  N0 = readRDS('../../../Vertebrate_Head/timepoints_rds/Amp_N0.RDS'),
  N2 = readRDS('../../../Vertebrate_Head/timepoints_rds/Amp_N2.RDS'),
  N5 = readRDS('../../../Vertebrate_Head/timepoints_rds/Amp_N5.RDS')
)

amp_list_renamed <- lapply(names(amp_list), function(x) {
  annot_sub<-coarse_annot[coarse_annot$stage==x,]
  x<-amp_list[[x]]
  x$celltype_fine<-as.character(Idents(x))
  x_r<-RenameIdents(x,setNames(annot_sub$coarse1,annot_sub$celltype_correct))
  x$celltype_coarse1<-as.character(Idents(x_r))
  x_r<-RenameIdents(x,setNames(annot_sub$coarse2,annot_sub$celltype_correct))
  x$celltype_coarse2<-as.character(Idents(x_r))
  return(x)
})

amp_list_sc <- lapply(amp_list_renamed, function(x) {
  #x$celltype <- as.character(Idents(x))
  x <- as.SingleCellExperiment(x)
  x <- x[rowSums(assay(x, 'counts')) > 0, colSums(assay(x, 'counts')) >
           0]
  x$Elav_plus <- assay(x, 'logcounts')['Elav', ] > 0
  x$celltype_fine_Eelav<-x$celltype_fine
  x$celltype_fine_Eelav[x$celltype_fine_Eelav == 'Ectoderm' & x$Elav_plus] <- 'Ectoderm_Elav'
  return(x)
})

names(amp_list_sc)<-names(amp_list)

detected_features <- Reduce(union, lapply(amp_list_sc, function(x) {
  x <- rownames(x)
  return(x)
}))

intersect_features <-
  Reduce(intersect, lapply(amp_list_sc, function(x) {
    x <- rownames(x)
    return(x)
  }))

max_cells <- max(sapply(amp_list_sc, ncol))

specific_mat <-
  matrix(rep_len(0, length(
    setdiff(detected_features, intersect_features)
  ) * max_cells), ncol = max_cells)
rownames(specific_mat) <-
  setdiff(detected_features, intersect_features)

amp_list_sc_full <- lapply(amp_list_sc, function(x) {
  new_sc_counts <-
    rbind(assay(x, 'counts'),
          specific_mat[setdiff(rownames(specific_mat), rownames(x)), 1:ncol(x)])
  new_sc_counts <- new_sc_counts[detected_features, ]
  new_sc <- SingleCellExperiment(assays = list(counts = new_sc_counts))
  colData(new_sc) <- colData(x)
  return(new_sc)
})

amp_list_sc_full <- lapply(names(amp_list_sc_full), function(x) {
  sce <- amp_list_sc_full[[x]]
  colnames(sce) <- paste0(colnames(sce), '_', x)
  colData(sce) <-
    colData(sce)[, c("nCount_RNA", "nFeature_RNA", "celltype_fine","celltype_coarse1","celltype_coarse2","celltype_fine_Eelav")]
  sce$stage <- x
  return(sce)
})

amp_sc_merged <- Reduce(cbind, amp_list_sc_full)

amp_sc_merged

amp_sc_merged$celltype_fine_stage<-paste0(amp_sc_merged$celltype_fine,'_',amp_sc_merged$stage)
amp_sc_merged$celltype_coarse1_stage<-paste0(amp_sc_merged$celltype_coarse1,'_',amp_sc_merged$stage)
amp_sc_merged$celltype_coarse2_stage<-paste0(amp_sc_merged$celltype_coarse2,'_',amp_sc_merged$stage)
amp_sc_merged$celltype_fine_Eelav_stage<-paste0(amp_sc_merged$celltype_fine_Eelav,'_',amp_sc_merged$stage)

writeH5AD(
  amp_sc_merged,
  file = 'amp_all.h5ad',
  colData = TRUE,
  rowData = FALSE,
  metadata = FALSE,
  version = "0.7.6"
)


system('chown rstudio:rstudio amp_all.h5ad')