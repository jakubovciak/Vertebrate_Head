library(Seurat)
library(zellkonverter)
library(Rhdf5lib)
library(HDF5Array)
library(scater)
library(dplyr)
library(scran)
#library(GenomicFeatures)


# prepare Amp data ####

coarse_annot<-read.csv('amp_coarse_annot.csv')

amp_rds<-list.files('../../../timepoints_rds',full.names = TRUE)

for (amp_path in amp_rds) {
  
  amp_name<-gsub(".RDS",'',basename(amp_path))
  stage<-gsub("Amp_","",amp_name)
  annot_sub<-coarse_annot[coarse_annot$stage==stage,]
  amp<-readRDS(amp_path)
  amp$celltype_fine<-as.character(Idents(amp))
  amp_r<-RenameIdents(amp,setNames(annot_sub$coarse1,annot_sub$celltype_correct))
  amp$celltype_coarse1<-as.character(Idents(amp_r))
  amp_r<-RenameIdents(amp,setNames(annot_sub$coarse2,annot_sub$celltype_correct))
  amp$celltype_coarse2<-as.character(Idents(amp_r))
  
  amp <- DietSeurat(
    amp,
    data = TRUE,
    scale.data = FALSE,
    dimreducs = NULL,
    graphs = NULL
  )
  amp<-as.SingleCellExperiment(amp)
  amp$Elav_plus <- assay(amp, 'logcounts')['Elav', ] > 0
  amp$celltype_fine_Eelav<-amp$celltype_fine
  amp$celltype_fine_Eelav[amp$celltype_fine_Eelav == 'Ectoderm' & amp$Elav_plus] <- 'Ectoderm_Elav'
  assay(amp,'logcounts')<-NULL
  

  
  amp$cluster_celltype<-paste(amp$seurat_clusters,amp$celltype,sep = '_')
  
  writeH5AD(
    amp,
    file = paste0('./',amp_name,'.h5ad'),
    colData = TRUE,
    rowData = FALSE,
    metadata = FALSE,
    version = "0.7.6"
  )
}

system('chown rstudio:rstudio *.h5ad')
