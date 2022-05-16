library(Seurat)
library(dplyr)
library(grid)
library(edgeR)
library(tibble)


setwd('./Working/Mouse_human_TG/20210401_mouse_TG_headache')


#load seurat object
seurat_mat=readRDS("Seurat.Rds")

dir = paste0(Sys.Date(),'_edgeR_IS_activated_naive_downsampledUMI')
dir.create(dir)
setwd(dir)

# get raw counts tables from Seurat object
cd = seurat_mat@assays$RNA@counts
dim(cd)

# get metadata
meta = seurat_mat@meta.data%>%rownames_to_column('V1')

# get cells from IS model that are activated by IEG scores
meta.activate=meta%>%filter(is.activated,model=='IS')

# count numbers
nCell.activate=meta.activate%>%
            group_by(class)%>%
            summarise(n=n())%>%
            as.data.frame()

# get naive cells that are not activated
meta.naive=meta%>%filter(!is.activated,model=='Naive')

print(nCell.activate)

# DE analysis
merge=NULL
# Repeat DE 10 times with randomly sampled naive populations
for(i in 1:10){
    print(i)

    meta.naive.select=NULL
    for(j in 1:nrow(nCell.activate)){
        #randomly select naive cells from the same seurat cluster as the activated cells
        meta.naive.select=rbind(meta.naive.select,meta.naive%>%filter(class==nCell.activate[j,'class'])%>%sample_n(nCell.activate[j,'n']))
    }
    
    condition_sub = rbind(meta.activate,meta.naive.select)%>%column_to_rownames('V1')

    cd_activate=cd[,meta.activate$V1]
    cd_naive=cd[,meta.naive.select$V1]

    # calculate factor for downsample
    nUMI.activate=meta.activate%>%pull(nCount_RNA)%>%mean()
    factor=(meta.naive.select%>%pull(nCount_RNA)%>%mean())/nUMI.activate
    
    # down sample counts table of activated cells to match the average UMI of naive cells
    cd_activate_down = apply(cd_activate,2, FUN=function(x) rbinom(length(x),x,factor))

    cd_sub=cbind(cd_naive,cd_activate_down)

    group = condition_sub$is.activated 
    group = as.factor(group)
    # IS activated cells are set as the baseline group, so DE compares naive cells to IS activated ones
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
    # reverse numbers. See line #66
    DEgenes$logFC = -DEgenes$logFC
    DEgenes$round=i
    print(dim(DEgenes))
    merge=rbind(merge,DEgenes%>%rownames_to_column('V1'))
}

data.table::fwrite(merge,file=paste0(Sys.Date(),'_edgeR_IS_activated_naive_downsampledUMI.csv'))





