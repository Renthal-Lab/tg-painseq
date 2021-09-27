library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(Hmisc)

setwd('/n/scratch3/users/l/ly116/Working/mTG_project_update')


min_genes=400
pc_num=30
res=1

date = format(Sys.time(),"%Y%m%d")
wd=paste(date,'Renthal_Sharma_Ryba_integrate',sep='_')
print(wd)
dir.create(wd)
setwd(wd)

# ============ Load data ============ 
# process NST
seurat_renthal=readRDS("../20201029_mouse_TG_level1_final/Seurat.Rds")
seurat_sharma=readRDS("../20201118_Sharma_TG_level1_final/Seurat.Rds")
seurat_ryba=readRDS("../20201118_Ryba_TG_level1_final/Seurat.Rds")

seurat_renthal$dataset='Renthal'
seurat_sharma$dataset='Sharma'
seurat_ryba$dataset='Ryba'


# ========= Step 2 =========
# Sctransform
seurat.list <- list(renthal=seurat_renthal,sharma=seurat_sharma,ryba=seurat_ryba)

# Integration of datasets
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:pc_num)
seurat_mat <- IntegrateData(anchorset = seurat.anchors, dims = 1:pc_num)

# ========= Step 3 =========
seurat_mat <- ScaleData(seurat_mat)

# PCA
seurat_mat <- RunPCA(seurat_mat, features = VariableFeatures(object = seurat_mat),verbose=F)

# dimension reduction and visualization
seurat_mat <- RunTSNE(seurat_mat, dims = 1:pc_num)
seurat_mat <- RunUMAP(seurat_mat, dims = 1:pc_num)

seurat_mat = AddMetaData(seurat_mat,Embeddings(seurat_mat[["tsne"]]),colnames(Embeddings(seurat_mat[["tsne"]])))
seurat_mat = AddMetaData(seurat_mat,Embeddings(seurat_mat[["umap"]]),colnames(Embeddings(seurat_mat[["umap"]])))

# ========= Step 4 =========
pdf("tSNE_dataset.pdf")
DimPlot(seurat_mat, reduction = "tsne",label = T, group.by='dataset')
dev.off()

pdf("UMAP_dataset.pdf")
DimPlot(seurat_mat, reduction = "umap",label = T, group.by='dataset')
dev.off()


markers =list(c('Rbfox3','Sparc','Rgs11','Meg3'),
        c('Atf3','Sox11','Sprr1a','Flrt3'),
        c('Fam19a4','Th'), #cLTMR
        c('Tac1',"Gpx3","Cartpt","Hpca","Trpm8","Calca","Trpv1","Scg2","Adcyap1"), #PEP
        c("Cd55","Mrgprd","Lpar3"), #NP
        c("Htr3a","Cplx2","Nptx1","Nefh","Hapln4","Pvalb","Cadps2","Ntrk2"), #NF
        c("Nppb","Sst","Il31ra"), #SST
        c('Apoe','Fabp7','Ednrb'),#Satglia
        c('Scn7a'),#Schwann_N
        c('Mpz','Mbp'),#Schwann_M,
        c('Igfbp7','Tinagl1','Rgs5','Myl9'),#vascular
        c('Gfap','Mlc1'), # astrocyte
        c('Mog','Hapln2'), # oligodendrocyte
        c('Dcn','Mgp','Pdgfra','Ngfr','Alpl'),#fibroblast
        c('Lyz2','Mrc1')# Immune
         )


# plot markers
dir.create('FeaturePlot')
index=0
for (i in markers)
{
  index=index+1
  i=i[i %in% rownames(seurat_mat)]
  if(length(i)==0){
    next()
  }
  print(i)
  pdf(paste0("FeaturePlot/FeaturePlot_tSNE_",index,".pdf"),heigh=ceiling(length(i)/2)*4,width=10)
  grid.draw(FeaturePlot(seurat_mat, i,reduction='tsne',ncol=2,coord.fixed=T))
  dev.off()

  pdf(paste0("FeaturePlot/FeaturePlot_UMAP_",index,".pdf"),heigh=ceiling(length(i)/2)*4,width=10)
  grid.draw(FeaturePlot(seurat_mat, i,reduction='umap',ncol=2,coord.fixed=T))
  dev.off()
}

# ========= Step 5 =========
# save files
print('saving files')
data.table::fwrite(seurat_mat@meta.data,file='meta.data.csv',row.names=TRUE)
saveRDS(seurat_mat, file = "Seurat.Rds")
print('done')


ident=seurat_mat@meta.data$subtype
names(ident)=rownames(seurat_mat@meta.data)
Idents(seurat_mat)= ident

pdf(paste0("tSNE_subtype.pdf"))
grid.draw(DimPlot(seurat_mat, reduction = "tsne",label = T,ncol=2))
dev.off()

pdf(paste0("UMAP_subtype.pdf"))
grid.draw(DimPlot(seurat_mat, reduction = "umap",label = T,ncol=2))
dev.off()

