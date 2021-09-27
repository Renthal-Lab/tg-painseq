library(Seurat)
library(dplyr)
library(grid)
library(edgeR)
library(tibble)


wd='/n/scratch3/users/l/ly116/Working/20210112_mouse_human_TG_project/20201029_human_TG_level1_final'

setwd(wd)

#load seurat obj
seurat_mat=readRDS("Seurat.Rds")

dir = paste0(Sys.Date(),'_edgeR_virusPos_virusNeg')
dir.create(dir)
setwd(dir)

cd = seurat_mat@assays$RNA@counts
dim(cd)

meta = seurat_mat@meta.data%>%rownames_to_column('V1')

counts.select=data.frame(t(as.matrix(cd[c('HSV1-LAT','GAPDH'),])))%>%
  rownames_to_column('V1')%>%
  select(V1,HSV1.LAT)

meta=left_join(meta,counts.select)%>%mutate(is_virus=HSV1.LAT>0)

meta.virusPos=meta%>%filter(is_virus)

nCell.virusPos=meta.virusPos%>%
            group_by(seurat_clusters)%>%
            summarise(n=n())%>%
            as.data.frame()

meta.virusNeg=meta%>%filter(!is_virus)

merge=NULL
for(i in 1:10){
    print(i)
    meta.virusNeg.select=NULL
    for(j in 1:nrow(nCell.virusPos)){
        meta.virusNeg.select=rbind(meta.virusNeg.select,meta.virusNeg%>%filter(seurat_clusters==nCell.virusPos[j,'seurat_clusters'])%>%sample_n(nCell.virusPos[j,'n']))
    }
    

    condition_sub = rbind(meta.virusPos,meta.virusNeg.select)%>%column_to_rownames('V1')

    cd_sub = cd[,rownames(condition_sub)]
    group = condition_sub$is_virus 
    group = as.factor(group)
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
    DEgenes$logFC = -DEgenes$logFC
    DEgenes$mean.1 = apply(cd_sub[rownames(DEgenes),condition_sub%>%rownames_to_column('V1')%>%filter(is_virus)%>%pull(V1)],1,mean,na.rm = TRUE)
    DEgenes$mean.2 = apply(cd_sub[rownames(DEgenes),condition_sub%>%rownames_to_column('V1')%>%filter(!is_virus)%>%pull(V1)],1,mean,na.rm = TRUE)
    DEgenes$pct.1 = apply(cd_sub[rownames(DEgenes),condition_sub%>%rownames_to_column('V1')%>%filter(is_virus)%>%pull(V1)]!=0,1,FUN = function(x) sum(x)/length(x))
    DEgenes$pct.2 = apply(cd_sub[rownames(DEgenes),condition_sub%>%rownames_to_column('V1')%>%filter(!is_virus)%>%pull(V1)]!=0,1,FUN = function(x) sum(x)/length(x))
    DEgenes$nCell.1 = nrow(condition_sub%>%filter(is_virus))
    DEgenes$nCell.2 = nrow(condition_sub%>%filter(!is_virus))
    DEgenes$round=i
    print(dim(DEgenes))
    merge=rbind(merge,DEgenes%>%rownames_to_column('V1'))
}


merge=merge%>%
            mutate(mean.1=round(mean.1,4),mean.2=round(mean.2,4),pct.1=round(pct.1,4),pct.2=round(pct.2,4))%>%
            select(V1,Log2FC=logFC,logCPM,F,PValue,mean.1,mean.2,pct.1,pct.2,nCell.1,nCell.2,round)

data.table::fwrite(merge,file=paste0(Sys.Date(),'_edgeR_virusPos_virusNeg','.csv'))





