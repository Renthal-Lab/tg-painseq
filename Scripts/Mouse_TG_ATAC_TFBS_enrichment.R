library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(tibble)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)


setwd('./Working/Mouse_human_TG/20210723_mTG_ATAC')


seurat_mat=readRDS('Seurat.Rds')

# ========== Add TFBS motif information =====================
# extract position frequency matrices for the motifs
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(all=T)
)

pfm=pfm[lapply(pfm, function(x) x@tags$species)=='Mus musculus']
# Homo sapiens Taxonomy ID: 9606
# Mus musculus Taxonomy ID: 10090

# add motif information
seurat_mat <- AddMotifs(
object = seurat_mat,
genome = BSgenome.Mmusculus.UCSC.mm10,
pfm = pfm
)

# ========== Calculate TFBS enrichment for each subtype =====================
# sub-type-specific peak results identifed by FindAllMarkers() in Seurat
da_peaks=data.table::fread('../20210723_mTG_ATAC_anchor_updated/subtypeMarkers.csv')


# format peaks
da_peaks=da_peaks%>%mutate(gene=stringr::str_replace(gene,':','-'))

# do it for each subtype
da_motifs=NULL
for (i in unique(da_peaks$subtype)){
	# get top differentially accessible peaks
	top.da.peak <- da_peaks%>%filter(subtype==i,p_val_adj < 0.05, avg_logFC >0.5)%>%pull(gene)
	print(i)
	print(length(top.da.peak))
	# calculate TFBS enrichment as the fold enrichment of counts of TFBS in sub-type-specific peaks compared to those in same number if randomly selected background peaks.
	# repeat for 10 times using different sets of background peaks
	for (n in 1:10){
		enriched.motifs <- FindMotifs(
		object = seurat_sub,
		features = top.da.peak,
		background = length(top.da.peak)
		)
		da_motifs=rbind(da_motifs,enriched.motifs%>%mutate(subtype=i,round=n,peaks=length(top.da.peak)))
	}
}

data.table::fwrite(da_motifs,file='motif_enrichment.csv')