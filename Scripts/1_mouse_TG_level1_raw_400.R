library(Seurat)
library(dplyr)
library(grid)
library(gridExtra)
library(tibble)
library(ggplot2)

setwd('/n/scratch3/users/l/ly116/Working/20210112_mouse_human_TG_project')

pc_num=20
min_genes=400
res=1.5

# make dirctory to save data
dir = paste0(format(Sys.time(),"%Y%m%d"),"_mouse_TG_level1_raw_",min_genes,"_",res)
dir.create(dir)
setwd(dir)

# read in counts table first
load('../counts/20191205_cd_s_tg_samples.Robj')

seurat_mat=CreateSeuratObject(counts = cd_s, project = "Mouse_TG", min.cells = 3, min.features = min_genes)
seurat_mat

# add sample meta data
file.rename('../sample.csv','sample.csv')
sample=data.table::fread('sample.csv')

meta=seurat_mat@meta.data%>%
    rownames_to_column('V1')%>%
    mutate(sample_name=stringr::str_remove(V1,'_bc.*'))%>%
    inner_join(sample)%>%
    column_to_rownames('V1')

dim(meta)
seurat_mat = AddMetaData(seurat_mat,meta)


# calculate mitohondrial genes for each cell
seurat_mat[["percent.mt"]] <- PercentageFeatureSet(seurat_mat, pattern = "^mt-")
seurat_mat <- subset(seurat_mat, subset = nFeature_RNA < 15000 & percent.mt < 5)

# plot QC
pdf(paste0("VinPlot_QC.pdf"),heigh=12,width=12)
VlnPlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Normalize data, normalization.method = "LogNormalize", scale.factor = 10000
seurat_mat <- NormalizeData(seurat_mat)

# ID variable genes, selection.method = "vst", nfeatures = 2000
seurat_mat <- FindVariableFeatures(seurat_mat)

# scale data
seurat_mat <- ScaleData(object = seurat_mat,vars.to.regress = c("nCount_RNA", "percent.mt"))

# PCA
seurat_mat <- RunPCA(seurat_mat, features = VariableFeatures(object = seurat_mat),verbose=F)

# dimension reduction and cluster
seurat_mat <- FindNeighbors(seurat_mat, dims = 1:pc_num)
seurat_mat <- FindClusters(seurat_mat, resolution = res)
print(levels(Idents(seurat_mat)))

# visualization
seurat_mat <- RunTSNE(seurat_mat, dims = 1:pc_num)
seurat_mat <- RunUMAP(seurat_mat, dims = 1:pc_num)

seurat_mat = AddMetaData(seurat_mat,Embeddings(seurat_mat[["tsne"]]),colnames(Embeddings(seurat_mat[["tsne"]])))
seurat_mat = AddMetaData(seurat_mat,Embeddings(seurat_mat[["umap"]]),colnames(Embeddings(seurat_mat[["umap"]])))

pdf("tSNE.pdf")
DimPlot(seurat_mat, reduction = "tsne",label = T)
dev.off()

pdf("UMAP.pdf")
DimPlot(seurat_mat, reduction = "umap",label = T)
dev.off()

pdf(paste0("tSNE_QC.pdf"),heigh=12,width=12)
FeaturePlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction='tsne',coord.fixed=T)
dev.off()

pdf(paste0("UMAP_QC.pdf"),heigh=12,width=12)
FeaturePlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction='umap',coord.fixed=T)
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
    pdf(paste0("FeaturePlot/FeaturePlot_tSNE_",index,".pdf"),heigh=ceiling(length(i)/2)*4,width=8)
    grid.draw(FeaturePlot(seurat_mat, i,reduction='tsne'))
    dev.off()

    pdf(paste0("FeaturePlot/FeaturePlot_UMAP_",index,".pdf"),heigh=ceiling(length(i)/2)*4,width=8)
    grid.draw(FeaturePlot(seurat_mat, i,reduction='umap'))
    dev.off()
}

# save files
print('saving files')
data.table::fwrite(seurat_mat@meta.data,file='meta.data.csv',row.names=TRUE)
saveRDS(seurat_mat, file = "Seurat.Rds")
print('done')


# save normalized counts tables
matrix=as.matrix(NormalizeData(seurat_mat,normalization.method='RC')@assays$RNA@data)
saveRDS(matrix,file='Scaled_10k_counts.Rds')

# find markers on a small subset to improve efficiency
cells.select=seurat_mat@meta.data%>%
      rownames_to_column('V1')%>%
      sample_frac(0.4)%>%
      pull(V1)


seurat_mat
seurat_mat=subset(seurat_mat,cells=cells.select)
seurat_mat

clusterMarkers <- FindAllMarkers(seurat_mat, only.pos = F, min.pct = 0.1, thresh.use = 0.5)
write.csv(clusterMarkers, file = "clusterMarkers.csv")
