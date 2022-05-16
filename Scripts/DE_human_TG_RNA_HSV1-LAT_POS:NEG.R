library(Seurat)
library(dplyr)
library(grid)
library(edgeR)
library(tibble)


setwd('./Working/Mouse_human_TG/20201029_human_TG_level1_final')


#load seurat object
seurat_mat=readRDS("Seurat.Rds")

dir = paste0(Sys.Date(),'_edgeR_HSV1_Pos_Neg')
dir.create(dir)
setwd(dir)

# get raw counts tables from Seurat object
cd = seurat_mat@assays$RNA@counts
dim(cd)

# get metadata
meta = seurat_mat@meta.data%>%rownames_to_column('V1')

# get HSV1_LAT positive cells
meta.virusPos=meta%>%filter(is_HSV1.LAT_pos)

# count numbers
nCell.virusPos=meta.virusPos%>%
            group_by(seurat_clusters)%>%
            summarise(n=n())%>%
            as.data.frame()

# get HSV1_LAT negative cells
meta.virusNeg=meta%>%filter(!is_virus)

# DE analysis
merge=NULL
# Repeat DE 10 times with randomly sampled negative populations
for(i in 1:10){
    print(i)
    meta.virusNeg.select=NULL
    for(j in 1:nrow(nCell.virusPos)){
        #randomly select negative cells from the same seurat cluster as the positive cells
        meta.virusNeg.select=rbind(meta.virusNeg.select,meta.virusNeg%>%filter(seurat_clusters==nCell.virusPos[j,'seurat_clusters'])%>%sample_n(nCell.virusPos[j,'n']))
    }
    

    condition_sub = rbind(meta.virusPos,meta.virusNeg.select)%>%column_to_rownames('V1')

    cd_sub = cd[,rownames(condition_sub)]
    group = condition_sub$is_virus 
    group = as.factor(group)
    # Positive cells are set as the baseline group, so DE compares negative cells to positive ones
    group = relevel(group, "TRUE")

    #edgeR
    y = DGEList(counts=cd_sub,group=group)
    y = calcNormFactors(y)
    design = model.matrix(~group, data=y$samples)
    y = estimateDisp(y,design)
    fit = glmQLFit(y,design)

    # quasi-likelihood (QL) F-test between naive vs later timepoint
    qlf = glmQLFTest(fit,coef=2)
    DEgenes = topTags(qlf, Inf)
    DEgenes = DEgenes$table
    # reverse numbers. See line #54
    DEgenes$logFC = -DEgenes$logFC
    DEgenes$round=i
    print(dim(DEgenes))
    merge=rbind(merge,DEgenes%>%rownames_to_column('V1'))
}

data.table::fwrite(merge,file=paste0(Sys.Date(),'_edgeR_HSV1_Pos_Neg.csv'))





