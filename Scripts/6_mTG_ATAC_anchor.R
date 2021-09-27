library(Seurat)
library(rtracklayer)
library(tibble)
library(dplyr)
library(gridExtra)
library(ggplot2)

setwd('/n/scratch3/users/l/ly116/Working/20210112_mouse_human_TG_project')

# Initiation
ATAC_species='Mouse'
RNA_species='Mouse'
pc_num=20

dir_reference=switch(ATAC_species,
                        'Human'="/n/scratch3/users/l/ly116/Reference_genome/RNA/GRCh38/genes/genes.gtf",
                        'Mouse'="/n/scratch3/users/l/ly116/Reference_genome/RNA/mm10/genes/genes.gtf")

# seq.levels: use chromosomes 1:22, X, and Y (human) or chromosomes 1:19, X, and Y (mouse).
seq_levels=switch(ATAC_species,
                        'Human'=paste0('chr',c(1:22, "X", "Y")),
                        'Mouse'=paste0('chr',c(1:19, "X", "Y")))


date = format(Sys.time(),"%Y%m%d")
wd=paste(date,'mTG_ATAC_anchor',sep='_')
print(wd)
dir.create(wd)
setwd(wd)
dir.create('ATAC_raw')
dir.create('ATAC_anchor')
dir.create('ATAC_filtered')


# =========== Step 1 ===========
# Read in peak matrix from 10X cellrange-ATAC output
dir_10X='../counts/20210518_mTG_ATAC_aggregate/outs'
peaks <- Read10X_h5(paste0(dir_10X,"/filtered_peak_bc_matrix.h5"))


# Create an imputed gene activity matrix from the peak matrix and GTF
# upstream: peaks that fall within gene bodies, or 2kb upstream of a gene, are considered
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = dir_reference, 
    seq.levels = seq_levels, upstream = 2000, verbose = TRUE, keep.sparse = TRUE)

# Convert human gene name to mouse gene name (e.g. ATF3 --> Atf3) for cross-species co-embedding
if(ATAC_species=='Human' && RNA_species=='Mouse'){
    rownames(activity.matrix)=capitalize(tolower(rownames(activity.matrix)))
}

# setup seurat object
seurat_atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10X_ATAC")
seurat_atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
seurat_atac$tech <- "atac"

# add 10X sequencing metrics to seurat object
meta <- read.table(paste0(dir_10X,"/singlecell.csv"), sep = ",", header = TRUE, row.names = 1, 
    stringsAsFactors = FALSE)
meta <- meta[colnames(seurat_atac), ]
seurat_atac <- AddMetaData(seurat_atac, metadata = meta)

#pre-process data
DefaultAssay(seurat_atac) <- "ACTIVITY"
seurat_atac <- FindVariableFeatures(seurat_atac)
seurat_atac <- NormalizeData(seurat_atac)
seurat_atac <- ScaleData(seurat_atac)

DefaultAssay(seurat_atac) <- "ATAC"
VariableFeatures(seurat_atac) <- names(which(Matrix::rowSums(seurat_atac) > 100))
seurat_atac <- RunLSI(seurat_atac, n = 50, scale.max = NULL)

# Visualization
seurat_atac <- RunTSNE(seurat_atac, reduction = "lsi", dims = 1:pc_num)
seurat_atac <- RunUMAP(seurat_atac, reduction = "lsi", dims = 1:pc_num)

seurat_atac = AddMetaData(seurat_atac,Embeddings(seurat_atac[["tsne"]]),colnames(Embeddings(seurat_atac[["tsne"]])))
seurat_atac = AddMetaData(seurat_atac,Embeddings(seurat_atac[["umap"]]),colnames(Embeddings(seurat_atac[["umap"]])))

# Save 
saveRDS(seurat_atac,file='ATAC_raw/Seurat.Rds')


# read in reference RNA seq data
seurat_rna=readRDS('../20201019_mouse_TG_level1_raw_400_1.5/Seurat.Rds')
seurat_rna$tech <- "RNA"

# read in subtype info
meta=data.table::fread('../20201029_mouse_TG_level1_final/meta.data.csv')
meta=meta%>%select(V1,class,cellID,subtype)%>%column_to_rownames('V1')

seurat_rna = AddMetaData(seurat_rna,meta,colnames(meta))
seurat_rna@meta.data[is.na(seurat_rna@meta.data)]='Doublet'

# downsample RNA data
nCell=nrow(seurat_atac@meta.data)*2
print(nCell)

cells.select=seurat_rna@meta.data%>%
        rownames_to_column('V1')%>%
        filter(model=='Naive')%>%
        sample_n(nCell)%>%
        pull(V1)
print(length(cells.select))

seurat_rna=subset(seurat_rna,cells=cells.select)
seurat_rna

# Identify anchors between the scATAC-seq dataset and the scRNA-seq dataset 
# and use these anchors to transfer the celltype labels we learned from the scRNA-seq data to the scATAC-seq cells.
transfer.anchors <- FindTransferAnchors(reference = seurat_rna, query = seurat_atac, features = VariableFeatures(object = seurat_rna), 
    reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = seurat_rna$subtype, 
    weight.reduction = seurat_atac[["lsi"]])
seurat_atac <- AddMetaData(seurat_atac, metadata = celltype.predictions)

# set cutoff for calling high confidence cells
cutoff=0.5
seurat_atac$hiConfi=ifelse(seurat_atac$prediction.score.max > cutoff,TRUE,FALSE)

table(seurat_atac$hiConfi)
seurat_atac@meta.data%>%group_by(predicted.id)%>%summarise(n.hiConfi=sum(hiConfi),pct.hiConfi=n.hiConfi/n())


# plot histograms of prediction id in each subtype
nCell_df=seurat_atac@meta.data%>%
    group_by(predicted.id)%>%
    summarise(nCell=n(),pct_hiConfi=round(sum(hiConfi)/nCell*100,2))

d1=ggplot(seurat_atac@meta.data,aes(x=prediction.score.max))+
    geom_histogram(bins=20)+
    geom_vline(xintercept=cutoff,color='red')+
    geom_text(data=nCell_df,aes(x=0.2,y=Inf,label=paste('Cell:',nCell)),vjust = 1)+
    geom_text(data=nCell_df,aes(x=0.5,y=Inf,label=paste('Pct_hc:',pct_hiConfi,'%')),vjust = 1)+
    facet_wrap(.~predicted.id)+
    theme_classic()

pdf('Histogram_prediction_id.pdf',height = 10,width = 13)
grid.arrange(d1)
dev.off()

# Save 
data.table::fwrite(seurat_atac@meta.data,file='ATAC_raw/meta.data.csv',row.names=TRUE)
saveRDS(seurat_atac,file='ATAC_raw/Seurat_atac.Rds')

# Co-embed
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(seurat_rna)
refdata <- GetAssayData(seurat_rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seurat_atac[["lsi"]])

# this line adds the imputed data matrix to the seurat_atac object
seurat_atac[["RNA"]] <- imputation
coembed <- merge(x = seurat_rna, y = seurat_atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)

coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)

# Visualization
coembed <- RunTSNE(coembed, dims = 1:pc_num)
coembed <- RunUMAP(coembed, dims = 1:pc_num)

coembed = AddMetaData(coembed,Embeddings(coembed[["tsne"]]),colnames(Embeddings(coembed[["tsne"]])))
coembed = AddMetaData(coembed,Embeddings(coembed[["umap"]]),colnames(Embeddings(coembed[["umap"]])))

d1 <- DimPlot(coembed, group.by = "tech")
d2 <- DimPlot(coembed, group.by = "subtype", label = TRUE, repel = TRUE)
d3 <- DimPlot(coembed, group.by = "hiConfi", label = TRUE, repel = TRUE, cols=c('TRUE'='green','FALSE'='red','NA'='transparent'))
d4 <- FeaturePlot(coembed, 'prediction.score.max',reduction = "umap",cols = c("grey","blue"))

pdf('UMAP_anchor.pdf',height = 10,width = 13)
grid.arrange(d1,d2,d3,d4,nrow=2)
dev.off()

# Save
data.table::fwrite(coembed@meta.data,file='ATAC_anchor/meta.data.csv',row.names=TRUE)
saveRDS(coembed,file='ATAC_anchor/Seurat.Rds')

# =========== Step 3 ===========
seurat_atac$subtype=seurat_atac$predicted.id
seurat_atac$prediciton.score=seurat_atac$prediction.score.max

# Subset high confidence cell and redo clustering and visualization
cells.select=seurat_atac@meta.data%>%filter(hiConfi,subtype != 'Doublet')

seurat_atac
seurat_atac=subset(seurat_atac,cells=rownames(cells.select))
seurat_atac


# clean up metadata by removing cell-type prediction scores
columns_remove=c(grep('prediction.score.',colnames(seurat_atac@meta.data),value=T),'hiConfi','cell_id','is_cell_barcode','predicted.id')

for( i in columns_remove){
    seurat_atac[[i]]=NULL
}


#pre-process data
DefaultAssay(seurat_atac) <- "ACTIVITY"
seurat_atac <- FindVariableFeatures(seurat_atac)
seurat_atac <- NormalizeData(seurat_atac)
seurat_atac <- ScaleData(seurat_atac)

DefaultAssay(seurat_atac) <- "ATAC"
VariableFeatures(seurat_atac) <- names(which(Matrix::rowSums(seurat_atac) > 100))
seurat_atac <- RunLSI(seurat_atac, n = 50, scale.max = NULL)


# Visualization
seurat_atac <- RunTSNE(seurat_atac, reduction = "lsi", dims = 1:pc_num)
seurat_atac <- RunUMAP(seurat_atac, reduction = "lsi", dims = 1:pc_num)

seurat_atac = AddMetaData(seurat_atac,Embeddings(seurat_atac[["tsne"]]),colnames(Embeddings(seurat_atac[["tsne"]])))
seurat_atac = AddMetaData(seurat_atac,Embeddings(seurat_atac[["umap"]]),colnames(Embeddings(seurat_atac[["umap"]])))

# set idents from subtype
Idents(object = seurat_atac) <- 'subtype'


pdf('tSNE_filtered.pdf',width=10,height=10)
grid.arrange(DimPlot(object = seurat_atac, label = TRUE,reduction = 'tsne') + NoLegend())
dev.off()


pdf('UMAP_filtered.pdf',width=10,height=10)
grid.arrange(DimPlot(object = seurat_atac, label = TRUE,reduction = 'umap') + NoLegend())
dev.off()

# save files
print('saving files')
data.table::fwrite(seurat_atac@meta.data,file='ATAC_filtered/meta.data.csv',row.names=TRUE)
saveRDS(seurat_atac, file = "ATAC_filtered/Seurat.Rds")
print('done')


clusterMarkers <- FindAllMarkers(seurat_mat)
write.csv(clusterMarkers, file = "clusterMarkers.csv")
