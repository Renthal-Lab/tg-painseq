# Human and mouse trigeminal ganglia cell atlas implicates multiple cell types in migraine 

## SUMMARY
The sensitization of trigeminal ganglion neurons contributes to primary headache disorders such as migraine, but the specific neuronal and non-neuronal trigeminal subtypes involved remain unclear. We thus developed a cell atlas in which human and mouse trigeminal ganglia are transcriptionally and epigenomically profiled at single-cell resolution. These data describe evolutionarily conserved and human-specific gene expression patterns within each trigeminal ganglion cell type, as well as the transcription factors and gene regulatory elements that contribute to cell-type-specific gene expression. We then leverage these data to identify trigeminal ganglion cell types that are implicated both by human genetic variation associated with migraine and two mouse models of headache. This trigeminal ganglion cell atlas improves our understanding of the cell types, genes, and epigenomic features involved in headache pathophysiology and establishes a rich resource of cell-type-specific molecular features to guide the development of more selective treatments for headache and facial pain.

## What's included in this directory
- [mouse_TG_RNA_cluster.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/mouse_TG_RNA_cluster.R) (Fig.1): Clustering and visualization of snRNA-seq data of mouse trigeminal ganglia.
- [Human_TG_RNA_cluster.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/Human_TG_RNA_cluster.R) (Fig.1): Clustering and visualization of snRNA-seq data of human trigeminal ganglia.
- [Human_mouse_TG_RNA_anchor.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/Human_mouse_TG_RNA_anchor.R) (Fig.2): Anchoring of snRNA-seq data of mouse trigeminal ganglia to human trigeminal ganglia and projection of human cell types to mouse nuclei.
- [DE_human_TG_RNA_HSV1-LAT_POS/NEG.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/DE_human_TG_RNA_HSV1-LAT_POS/NEG.R) (Fig.4): Differential expression analysis of snRNA-seq data of human trigeminal ganglia comparing HSV1-LAT+ nuclei to HSV1-LAT- nuclei using edgeR.
- [Mouse_TG_ATAC_cluster.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/Mouse_TG_ATAC_cluster.R) (Fig.5): Clustering and visualization of snATAC-seq data of mouse trigeminal ganglia by multi-modal anchoring to snRNA-seq data of mouse trigeminal ganglia.
- [Mouse_TG_ATAC_TFBS_enrichment.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/Mouse_TG_ATAC_TFBS_enrichment.R) (Fig.5): Enrichment analysis of transcription factor binding sites on peaks of human snATAC-seq cell types.
- [Mouse_TG_ATAC_peak_RNA_gene_correlation.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/Mouse_TG_ATAC_peak_RNA_gene_correlation.R) (Fig.6): Correlation analysis of the chromatin accessibility of cell-type-specific peaks in snATAC-seq with expression of cell-type-specific genes in snRNA-seq of mouse cell types.
- [Mouse_TG_ATAC_migraine_SNP_enrichment.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/Mouse_TG_ATAC_migraine_SNP_enrichment.R) (Fig.6): Enrichment of migraine-associated single nucleotide polymorphism in cell-type-specific peaks of snATAC-seq in mouse cell tpes .
- [Mouse_TG_RNA_IEG_score.Rmd](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/Mouse_TG_RNA_IEG_score.Rmd) (Fig.7): Identification of activated nuclei in mouse models of headache based on expression of immediate early genes.
- [DE_mouse_TG_RNA_headache_naive.R](https://github.com/Renthal-Lab/tg-painseq/blob/main/Scripts/DE_mouse_TG_RNA_headache_naive.R) (Fig.7): Differential expression analysis of snRNA-seq data of mouse trigeminal ganglia comparing activated nuclei from headache models to nuclei of naive animals.

## Data availability
-   [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197289)
-   [Trigeminal Ganglion Cell Atlas] (https://painseq.shinyapps.io/tg-painseq/)


## Cite our paper

[Yang et al., Human and mouse trigeminal ganglia cell atlas implicates multiple cell types in migraine, Neuron (2022)](https://doi.org/10.1016/j.neuron.2022.03.003)

## Related packages
-   [Seurat](https://github.com/satijalab/seurat)
-   [Signac](https://github.com/timoast/signac)
-   [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
