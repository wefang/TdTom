# Original analysis conducted with Seurat 2.0
# The clustering is used in the new analysis in Seurat 3.0
# TSNE plots etc. are replaced with UMAP and Vairable Gene Selection etc are carried out again with newer versions
d77=read.table("OPC_d77_dge.txt.gz", header = T, row.names=1)
d89=read.table("OPC_d89_dge.txt.gz", header = T, row.names=1)
d104=read.table("OPC_d104_dge.txt.gz", header = T, row.names=1)
d77=CreateSeuratObject(d77, min.cells=1, min.genes=1)
d89=CreateSeuratObject(d89, min.cells=1, min.genes=1)
d104=CreateSeuratObject(d104, min.cells=1, min.genes=1)
OPCs_tdTom=MergeSeurat(d77, d89, min.cells=1, min.genes=1, do.normalize=F, add.cell.id1="d77", add.cell.id2="d89")
OPCs_tdTom=MergeSeurat(OPCs_tdTom, d104, min.cells=3, min.genes=250, do.normalize=T, add.cell.id2="d104")

OPCs_tdTom=StashIdent(OPCs_tdTom, save.name="age")
VlnPlot(OPCs_tdTom, c("nUMI", "nGene"))
mito.genes=grep("^MT-", rownames(OPCs_tdTom@data), value=T)
percent.mito=Matrix::colSums(OPCs_tdTom@raw.data[mito.genes,])/Matrix::colSums(OPCs_tdTom@raw.data)
OPCs_tdTom=AddMetaData(OPCs_tdTom, metadata=percent.mito, col.name="percent.mito")
GenePlot(OPCs_tdTom, "percent.mito", "nUMI")
OPCs_tdTom=FilterCells(OPCs_tdTom, subset.names=c("nUMI", "percent.mito"), low.thresholds=c(350, -Inf),
high.thresholds=c(30000, 0.2))

# OPCs_tdTom=FindVariableGenes(OPCs_tdTom, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff=0.0125,
# x.high.cutoff=3, y.cutoff=0.5)
# OPCs_tdTom=FindVariableGenes(OPCs_tdTom, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff=0.025,
# x.high.cutoff=3, y.cutoff=0.5)
# OPCs_tdTom=FindVariableGenes(OPCs_tdTom, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff=0.025, x.high.cutoff=3, y.cutoff=1)
OPCs_tdTom=FindVariableGenes(OPCs_tdTom, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff=0.025, x.high.cutoff=3, y.cutoff=0.75)
OPCs_tdTom=RunPCA(OPCs_tdTom, pc.genes=OPCs_tdTom@var.genes, pcs.print=1:5, genes.print=5)
OPCs_tdTom=ScaleData(OPCs_tdTom)
OPCs_tdTom=ProjectPCA(OPCs_tdTom)
PCHeatmap(OPCs_tdTom, pc.use=1:12, cells.use=500, do.balanced=T)
PCHeatmap(OPCs_tdTom, pc.use=1:6, cells.use=500, do.balanced=T) 
PCHeatmap(OPCs_tdTom, pc.use=7:12, cells.use=500, do.balanced=T)
JackStrawPlot(OPCs_tdTom, PCs=1:20)

OPCs_tdTom=FindClusters(OPCs_tdTom, dims.use=1:16, resolution=1.2, reduction.type="pca", print.output=0, save.SNN=T)
OPCs_tdTom=RunTSNE(OPCs_tdTom, dims.use=1:16)
TSNEPlot(OPCs_tdTom, do.label=T)