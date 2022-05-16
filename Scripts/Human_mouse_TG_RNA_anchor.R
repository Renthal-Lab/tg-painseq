library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(Hmisc)

setwd('./Working/Mouse_human_TG')

# In this script, we are going to anchor human TG snRNA-seq data to mouse TG snRNA-seq data and transfer the mouse cell type labels to each of the human cells
# Adapted from Seurat tutorial: https://satijalab.org/seurat/articles/integration_mapping.html

# ============ Load data ============ 
seurat_human=readRDS("20201029_human_TG_level1_final/Seurat.Rds") # query
seurat_mouse=readRDS("20201029_mouse_TG_level1_final/Seurat.Rds") # reference


date = format(Sys.time(),"%Y%m%d")
wd=paste(date,'hTG_mTG_anchor',sep='_')
print(wd)
dir.create(wd)
setwd(wd)


# Names of human genes are converted by captalizing gene names
matrix_human=seurat_human@assays$RNA@counts
rownames(matrix_human)=capitalize(tolower((rownames(matrix_human))))

meta_human=seurat_human@meta.data

# create a seurat object from human dataset that is compatible with mouse dataset
seurat_human <- CreateSeuratObject(counts = matrix_human, project = "Human")
seurat_human

#preprocessing of human data
seurat_human <- AddMetaData(seurat_human, metadata = meta_human)
seurat_human <- NormalizeData(seurat_human)
seurat_human <- FindVariableFeatures(seurat_human)
seurat_human <- ScaleData(object = seurat_human, vars.to.regress = c("nCount_RNA", "percent.mt"))
seurat_human <- RunPCA(seurat_human, features = VariableFeatures(object = seurat_human),verbose=F)


# ============ projection of reference data onto query object ============ 
# We use all default parameters here for identifying anchors
transfer.anchors <- FindTransferAnchors(reference = seurat_mouse, query = seurat_human, features = VariableFeatures(object = seurat_mouse), 
    reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

# TransferData returns a matrix with predicted IDs and prediction scores, which are written into meta.data
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = seurat_mouse$subtype, 
    weight.reduction = seurat_human[["pca"]])
seurat_human <- AddMetaData(seurat_human, metadata = celltype.predictions)


# Cells with low anchoring scores (<=0.5) are either low-quality cells or cells with ambigous identity.
# They are marked as low-confidence cells and removed later.
cutoff=0.5
seurat_human$hiConfi=ifelse(seurat_human$prediction.score.max > cutoff,'TRUE','FALSE')

# summary report of anchoring results
table(seurat_human$predicted.id == seurat_human$subtype)
seurat_human@meta.data%>%mutate(agree.subtype=(predicted.id==subtype))%>%group_by(subtype)%>%summarise(n=sum(agree.subtype)/n())


# ============ Intergration of two datasets for visualization ============ 

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(seurat_mouse)
refdata <- GetAssayData(seurat_mouse, assay = "RNA", slot = "data")[genes.use, ]

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

# summary plots
d1 <- DimPlot(coembed, group.by = "species")
d2 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Mouse')%>%pull(V1),group.by = "subtype", label = TRUE, repel = TRUE)+labs(title='Mouse subtype')

d3 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Human')%>%pull(V1),group.by = "subtype", label = TRUE, repel = TRUE)+labs(title='Human subtype')
d4 <- DimPlot(coembed, cells= coembed@meta.data%>%rownames_to_column('V1')%>%filter(species=='Human')%>%pull(V1),group.by = "predicted.id", label = TRUE, repel = TRUE)+labs(title='Human predicted.id')

d5 = FeaturePlot(coembed, 'prediction.score.max',reduction = "umap",cols = c("grey","blue"))
d6 <- DimPlot(coembed, group.by = "hiConfi", label = TRUE, repel = TRUE, cols=c('TRUE'='green','FALSE'='red','NA'='transparent'))

pdf('UMAP_coembed.pdf',height = 15,width = 13)
grid.arrange(d1,d2,d3,d4,d5,d6,ncol=2)
dev.off()

saveRDS(coembed,file='Seurat_coembed.Rds')
