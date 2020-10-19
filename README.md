# TdTom
* Single cell RNA-seq analysis for human stem cell-derived oligodendrocyte lineage cells (GSE146373)

- Older analysis conducted for early versions of the paper.
	1. **seurat2_old_analysis.R** original analysis conducted with Seurat 2.0.
	2. **monocle_old_analysis.R** pseudo-temporal trajectory analysis conducted with monocle 2.0. (Fig 7-8)
- Final analysis conducted for the final version of the paper.
	1. **reanalyze_tdTom.R** standard Seurat 3.1.5 pipeline is used. Because the Louvain clustering produced similar results as the previous TSNE clustering, TSNE clusters are displayed on the UMAP plot. Marker gene lists are used for preranked GSEA and 
- Processing for data integration 
	1. **jakel.R** [data source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118257)
	2. **schimer19.R** [data source](https://cells.ucsc.edu/?ds=ms)
	3. **zhong.R** [data source](https://db.cngb.org/HCL/gallery.html?tissue=Fetal-Brain_Zhong)
- Integration and label transfer
	1. **integrate_multiple.R** integration of the three datasets above with data in the paper. (Fig 4)
- Additional supplementary data analysis:
	1. **hESC_10x.R** analysis for independently differentiated cells from day 85 (Fig S7)
	2. **ipsc.R** analysis for day 60 iPSC cells (Fig S8)