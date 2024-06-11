library(URD)
library(zellkonverter)
library(Rhdf5lib)
library(HDF5Array)
library(scater)
library(dplyr)
library(scran)

# prepare ZF data ####

zf_urd<-readRDS('URD_Zebrafish_Object.rds')

zf_segment<-read.table('URD_Dropseq_Meta.txt',row.names = 1,header = TRUE) %>% .[-1,1:2]

zf_lineage<-read.table('URD_Dropseq_Meta.txt',row.names = 1,header = TRUE) %>% 
  .[-1,] %>% 
  mutate(across(starts_with("Lineage"),as.logical)) %>% 
  select(-c(1:2)) %>% 
  mutate(celltype=apply(.,1,function(l_row){
    if(sum(l_row)==1){
      return(colnames(.)[l_row])
    } else {
      return("unassigned")
    }
  }))

zf_sce<-SingleCellExperiment(assays=list(counts=zf_urd@count.data))
colData(zf_sce)<-DataFrame(cbind(zf_urd@meta,celltype=as.character(zf_lineage$celltype)))
zf_sce$segment<-zf_segment$Segment
zf_sce$segment_celltype<-paste(zf_sce$segment,zf_sce$celltype,sep = '_')
zf_sce$stage_fine<-paste0(zf_sce$HPF,"_",zf_sce$STAGE)
zf_sce$segment_celltype_stage<-paste0(zf_sce$segment_celltype,"_",zf_sce$stage_fine)
zf_sce$celltype_stage<-paste0(zf_sce$celltype,"_",zf_sce$stage_fine)

zf_sce$celltype_coarse


rownames(zf_sce)[rownames(zf_sce)=='DEDD1']<-'DEDD1_1'

rownames_orig<-rownames(zf_sce)
saveRDS(rownames_orig,'rownames_orig.RDS')

rownames(zf_sce)<-tolower(rownames(zf_sce))


zf_sce_50to12<-zf_sce[,zf_sce$stage_fine%in%c("5.3_ZF50","7_ZF60","8_ZF75","9_ZF90","10_ZFB","11_ZF3S","12_ZF6S")]
zf_sce_50to12<-zf_sce_50to12[rowSums(assay(zf_sce_50to12,'counts'))>0,]
#zf_sub$stage_fine<-droplevels(zf_sub$stage_fine)

writeH5AD(
  zf_sce_50to12,
  file = 'zf_farrell_50to12.h5ad',
  colData = TRUE,
  rowData = FALSE,
  metadata = FALSE,
  version = "0.7.6"
)


for (stagehpf in unique(zf_sce$stage_fine)) {
  zf_sub<-zf_sce[,zf_sce$stage_fine%in%stagehpf]
  zf_sub<-zf_sub[rowSums(assay(zf_sub,'counts'))>0,]
  #zf_sub$stage_fine<-droplevels(zf_sub$stage_fine)
  
  writeH5AD(
    zf_sub,
    file = paste0('zf_farrell',stagehpf,'.h5ad'),
    colData = TRUE,
    rowData = FALSE,
    metadata = FALSE,
    version = "0.7.6"
  )

}

