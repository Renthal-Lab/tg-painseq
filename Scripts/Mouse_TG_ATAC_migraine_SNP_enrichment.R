library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(grid)
library(tibble)
library(GenomicRanges)
library(reshape)


setwd('./Working/Mouse_human_TG/20210526_mTG_ATAC')


seurat_mat=readRDS('Seurat.Rds')

# read in lifeOver SNPs, extend 500bp each side from center
ranges=data.table::fread('SNP_liftover_hg19_mm10_ext_5bp_eachside.bed',header=F)
ranges=ranges[[1]]
ranges=StringToGRanges(ranges, sep = c(":", "-"))
ranges=Extend(ranges, upstream = 500, downstream = 500, from.midpoint = T)

# calculate fraction of counts in each SNP region in each cell using FractionCountsInRegion() in Signac
# Reference of FractionCountsInRegion(): https://satijalab.org/signac/reference/fractioncountsinregion
counts=NULL
for(i in 1:length(ranges)){
  counts_temp=FractionCountsInRegion(seurat_mat,ranges[c(i,i)])
  counts=rbind(counts,counts_temp)
}

# format results
rownames(counts)=GRangesToString(ranges,sep = c(":", "-"))

# melt the data table, where row is cell ids, column is SNP coordinates, and the value in table is fraction of counts
# into a data frame of three columns: SNP coordinates, cell ids, fraction of counts
counts_df=melt(counts)
colnames(counts_df)=c('SNP','V1','value')

# summarise fraction of counts by cell type
data=seurat_mat@meta.data%>%
      rownames_to_column('V1')%>%
      select(V1,subtype)

data=data%>%left_join(counts_df)

data.table::fwrite(data,file='FractionCountsInSNP_subtype.csv')

