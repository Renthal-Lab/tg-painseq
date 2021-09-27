library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)

setwd('/n/scratch3/users/l/ly116/Working/20210112_mouse_human_TG_project')

# ============ Load data ============ 
seurat_human=readRDS("20201029_human_TG_level1_final/Seurat.Rds") # query
seurat_mouse=readRDS("20201029_mouse_TG_level1_final/Seurat.Rds") # reference

seurat_human
seurat_mouse

date = format(Sys.time(),"%Y%m%d")
wd=paste(date,'hTG_mTG_anchor',sep='_')
print(wd)
dir.create(wd)
setwd(wd)


# create a seurat objec from human dataset that is compatible with mouse dataset
# Gene names are converted by captalizing human gene names
matrix_human=seurat_human@assays$RNA@counts
rownames(matrix_human)=capitalize(tolower((rownames(matrix_human))))

meta_human=seurat_human@meta.data

seurat_human <- CreateSeuratObject(counts = matrix_human, project = "Human")
seurat_human

seurat_human <- AddMetaData(seurat_human, metadata = meta_human)
seurat_human <- NormalizeData(seurat_human)
seurat_human <- FindVariableFeatures(seurat_human)
seurat_human <- ScaleData(object = seurat_human, vars.to.regress = c("nCount_RNA", "percent.mt"))
seurat_human <- RunPCA(seurat_human, features = VariableFeatures(object = seurat_human),verbose=F)

# subset mouse dataset to ensure same subtypes with human
cells.select=seurat_mouse@meta.data%>%
			rownames_to_column('V1')%>%
			filter(subtype %in% unique(meta_human$subtype))%>%
			pull(V1)

seurat_mouse
seurat_mouse <- subset(seurat_mouse, cells=cells.select)
seurat_mouse


# ============ Co-embedding subtype ============ 
transfer.anchors <- FindTransferAnchors(reference = seurat_mouse, query = seurat_human, features = VariableFeatures(object = seurat_mouse), 
    reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = seurat_mouse$subtype, 
    weight.reduction = seurat_human[["pca"]])
seurat_human <- AddMetaData(seurat_human, metadata = celltype.predictions)
head(seurat_human@meta.data)

seurat_human$hiConfi=ifelse(seurat_human$prediction.score.max > cutoff,'TRUE','FALSE')
table(seurat_human$predicted.id == seurat_human$subtype)
seurat_human@meta.data%>%mutate(agree.subtype=(predicted.id==subtype))%>%group_by(subtype)%>%summarise(n=sum(agree.subtype)/n())

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(seurat_mouse)
refdata <- GetAssayData(seurat_mouse, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seurat_human[["pca"]])

# this line adds the imputed data matrix to the seurat_human object
seurat_human[["RNA"]] <- imputation
coembed <- merge(x = seurat_mouse, y = seurat_human)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunTSNE(coembed, dims = 1:20)
coembed <- RunUMAP(coembed, dims = 1:20)

coembed = AddMetaData(coembed,Embeddings(coembed[["tsne"]]),colnames(Embeddings(coembed[["tsne"]])))
coembed = AddMetaData(coembed,Embeddings(coembed[["umap"]]),colnames(Embeddings(coembed[["umap"]])))

coembed$species=ifelse(!is.na(coembed$predicted.id), 'Human', 'Mouse')


d1 <- DimPlot(coembed, group.by = "species")
d2 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Mouse')%>%pull(V1),group.by = "subtype", label = TRUE, repel = TRUE)+labs(title='Mouse subtype')

d3 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Human')%>%pull(V1),group.by = "subtype", label = TRUE, repel = TRUE)+labs(title='Human subtype')
d4 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Human')%>%pull(V1),group.by = "predicted.id", label = TRUE, repel = TRUE)+labs(title='Human predicted.id')

d5 = FeaturePlot(coembed, 'prediction.score.max',reduction = "umap",cols = c("grey","blue"))
d6 <- DimPlot(coembed, group.by = "hiConfi", label = TRUE, repel = TRUE, cols=c('TRUE'='green','FALSE'='red','NA'='transparent'))

pdf('UMAP_coembed_subtype.pdf',height = 15,width = 13)
grid.arrange(d1,d2,d3,d4,d5,d6,ncol=2)
dev.off()

data.table::fwrite(coembed@meta.data,file='meta.data.csv',row.names=TRUE)
saveRDS(coembed,file='Seurat_subtype.Rds')



# ============ Co-embedding subtype ============ 
transfer.anchors <- FindTransferAnchors(reference = seurat_mouse, query = seurat_human, features = VariableFeatures(object = seurat_mouse), 
    reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = seurat_mouse$cellID, 
    weight.reduction = seurat_human[["pca"]])
seurat_human <- AddMetaData(seurat_human, metadata = celltype.predictions)
head(seurat_human@meta.data)

seurat_human$hiConfi=ifelse(seurat_human$prediction.score.max > cutoff,'TRUE','FALSE')
table(seurat_human$predicted.id == seurat_human$cellID)
seurat_human@meta.data%>%mutate(agree.cellID=(predicted.id==cellID))%>%group_by(cellID)%>%summarise(n=sum(agree.cellID)/n())

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(seurat_mouse)
refdata <- GetAssayData(seurat_mouse, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seurat_human[["pca"]])

# this line adds the imputed data matrix to the seurat_human object
seurat_human[["RNA"]] <- imputation
coembed <- merge(x = seurat_mouse, y = seurat_human)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunTSNE(coembed, dims = 1:20)
coembed <- RunUMAP(coembed, dims = 1:20)

coembed = AddMetaData(coembed,Embeddings(coembed[["tsne"]]),colnames(Embeddings(coembed[["tsne"]])))
coembed = AddMetaData(coembed,Embeddings(coembed[["umap"]]),colnames(Embeddings(coembed[["umap"]])))

coembed$species=ifelse(!is.na(coembed$predicted.id), 'Human', 'Mouse')


d1 <- DimPlot(coembed, group.by = "species")
d2 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Mouse')%>%pull(V1),group.by = "cellID", label = TRUE, repel = TRUE)+labs(title='Mouse cellID')

d3 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Human')%>%pull(V1),group.by = "cellID", label = TRUE, repel = TRUE)+labs(title='Human cellID')
d4 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Human')%>%pull(V1),group.by = "predicted.id", label = TRUE, repel = TRUE)+labs(title='Human predicted.id')

d5 = FeaturePlot(coembed, 'prediction.score.max',reduction = "umap",cols = c("grey","blue"))
d6 <- DimPlot(coembed, group.by = "hiConfi", label = TRUE, repel = TRUE, cols=c('TRUE'='green','FALSE'='red','NA'='transparent'))

pdf('UMAP_coembed_cellID.pdf',height = 15,width = 13)
grid.arrange(d1,d2,d3,d4,d5,d6,ncol=2)
dev.off()

data.table::fwrite(coembed@meta.data,file='meta.data.csv',row.names=TRUE)
saveRDS(coembed,file='Seurat_subtype.Rds')