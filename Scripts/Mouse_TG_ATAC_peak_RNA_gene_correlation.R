library(dplyr)
library(tibble)
library(biomaRt)
library(GenomicRanges)
library(Signac)

setwd('./Working/Mouse_human_TG')

# In this script, we will perform correlation analysis of the chromatin accessibility of cell-type-specific peaks in snATAC-seq with expression of cell-type-specific genes in snRNA-seq of mouse cell types.
# To do that, we will:
# 1. Calculate the average expression of cell-type-specific genes in each of the cell types from mouse TG snRNA-seq data
# 2. Calculate the average chromatin accessibility of cell-type-specific peaks in each of the cell types from mouse TG snATAC-seq data
# 3. Calculate the pair-wise correlation between any cell-type-specific genes to any cell-type-specific peaks
# 4. Get the genomic positions of each of the cell-type-specific genes and cell-type-specific peaks


# 1. Calculate the average expression of cell-type-specific genes in each of the cell types from mouse TG snRNA-seq data
seurat_RNA=readRDS('/20201029_mouse_TG_level1_final/Seurat.rds') # Seurat object for snRNA-seq data
gene_DE_df=data.table::fread('/20201029_mouse_TG_level1_final/subytpeMarkers.csv') # Cell-type-specific gene list by FindAllMarkers()

# select cell-type-specific genes
gene_DE_sig=gene_DE_df%>%filter(p_val_adj<0.05,avg_logFC>0.5)%>%pull(gene)

# calculate average expression per cell type, and filter for cell-type-specific genes
gene=AverageExpression(seurat_RNA)$RNA%>%column_to_rownames("V1")%>%filter(V1 %in% gene_DE_sig)

# 2. Calculate the average chromatin accessibility of cell-type-specific peaks in each of the cell types from mouse TG snATAC-seq data
seurat_ATAC=readRDS('/20210202_mouse_TG_ATAC_final/Seurat.rds') # Seurat object for snRNA-seq data
peak_DE_df=data.table::fread('/20210202_mouse_TG_ATAC_final/subytpeMarkers.csv') # Cell-type-specific gene list by FindAllMarkers()

# select cell-type-specific peaks
peak_DE_sig=peak_DE_df%>%filter(p_val_adj<0.05,avg_logFC>0.5)%>%pull(peak)

# calculate average chromatin acessibility per cell type, and filter for cell-type-specific genes
peak=AverageExpression(seurat_ATAC)$ATAC%>%column_to_rownames("V1")%>%filter(V1 %in% peak_DE_sig)

# 3. Calculate the pair-wise correlation between any cell-type-specific genes to any cell-type-specific peaks
# select common cell types in both RNA-seq and ATAC-seq data
celltypes=intersect(colnames(peak),colnames(gene))

dim(peak)
dim(gene)
# select common cell types
peak=peak[,celltypes]
gene=gene[,celltypes]

# filter out features with low level of expression/accessibility
peak=t(peak[rowMeans(peak)>0.05,])
gene=t(gene[rowMeans(gene)>0.05,])

dim(peak)
dim(gene)

correlation=round(cor(peak,gene),3) # round up results to save spaces
saveRDS(correlation,'correlation.Rds')

# melt the data table, where row is peak, column is gene, and the value in table is correlation coeffcients,
# into a data frame of three columns: peak, gene, fcorrelation coeffcients
correlation.melt=as.data.frame(correlation)%>%
		rownames_to_column('peak')%>%
		melt(id.vars='peak',measure.vars=colnames(correlation))%>%
		dplyr::rename(gene=variable,cor=value)


# 4. Get the genomic positions of each of the cell-type-specific genes and cell-type-specific peaks
# A function to get genomic ranges for gene name input
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]]),
    "gene_name"
  ]
  colnames(x = collapsed) <- c(
    "gene_name", "seqnames", "start", "end", "strand", "gene_biotype"
  )
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

# get genomic information for genes
gene_coor=CollapseToLongestTranscript(Annotation(seurat_RNA))

# format data
gene_coor=gene_coor%>%
		as_tibble()%>%
		select(chromosome=seqnames,gene.start=start,gene.end=end,gene.strand=strand,gene=gene_name)%>%
		filter(gene %in% colnames(correlation))%>%
		mutate(gene.center=(gene.start+gene.end)/2)%>%
		select(gene,chromosome,gene.start,gene.end,gene.strand,gene.center)


# get genomic information for peaks
peak=stringr::str_split(rownames(correlation),'-',simplify=T)%>%
		as_tibble()%>%
		dplyr::rename(chromosome=V1,peak.start=V2,peak.end=V3)%>%
		mutate(peak=paste(chromosome,peak.start,peak.end,sep='-'),peak.start=as.numeric(peak.start),peak.end=as.numeric(peak.end))%>%
		mutate(peak.center=(peak.start+peak.end)/2)%>%
		select(peak,chromosome,peak.start,peak.end,peak.center)

# join data
correlation.merge=correlation.melt%>%
	inner_join(peak)%>%
	inner_join(gene) # join by gene and chromosome, so that only peak-gene pair on the same chromosome will be kept

dim(correlation.merge)

# calculate distance
correlation.merge$distance=ifelse(correlation.merge$gene.strand=='+',
				correlation.merge$peak.center-correlation.merge$gene.center,
				correlation.merge$gene.center-correlation.merge$peak.center)


fwrite(correlation.merge,file='CTSpeak_CTSgene_correlation_distance.csv')

