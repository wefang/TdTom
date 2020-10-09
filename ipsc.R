library(Seurat)
ipsc = CreateSeuratObject(Read10X("extra_data/gsm4672507/ipsc/"), min.cells = 10, min.features = 100)
mito.genes <- grep(pattern = "^MT-", x = rownames(ipsc), value = TRUE)
percent.mito <- Matrix::colSums(ipsc@assays$RNA@counts[mito.genes, ])/Matrix::colSums(ipsc@assays$RNA@counts)
ipsc <- AddMetaData(ipsc, metadata = percent.mito, col.name = "percent.mito")
ipsc = ipsc[, ipsc$percent.mito < 0.09]
VlnPlot(ipsc, c("nCount_RNA", "nFeature_RNA", "percent.mito"))

ipsc = FindVariableFeatures(ipsc, nfeatures = 2000)
ipsc = NormalizeData(ipsc)
ipsc = ScaleData(ipsc)
ipsc = RunPCA(ipsc)
ipsc = RunUMAP(ipsc, dims = 1:10)
ipsc = FindNeighbors(ipsc)
ipsc = FindClusters(ipsc)
cl_order = c(2, 6, 0 ,3 ,1, 4, 5)
cl_labels = c("CyP1", "CyP2", "OPC1", "OPC2", "OPC3", "OPC4", "OPC5")
cl_col = c(brewer.pal(4, "Purples")[3:4],
           brewer.pal(6, "RdPu")[2:6])
Idents(ipsc) <- factor(x = ipsc$seurat_clusters, levels = cl_order, labels = cl_labels)
UMAPPlot(ipsc, label = T, label.size = 10, pt.size = 3.0, repel = T) + NoLegend() + scale_color_manual(values = cl_col)
ipsc$source = "iPSC-OPCs"

g1s_genes = readLines("g1s_genes.txt")
g2m_genes = readLines("g2m_genes.txt")
ipsc = CellCycleScoring(ipsc, s.features = g1s_genes, g2m.features = g2m_genes)
FeaturePlot(ipsc, features = "S.Score") + FeaturePlot(ipsc, features = "G2M.Score")
ggsave("plot_extra/ipsc/psc_cell_cylce.pdf")

ipsc_markers = FindAllMarkers(ipsc)
readr::write_csv(ipsc_markers, "plot_extra/ipsc_markers.csv")
make_rnk(ipsc_markers, path = "plot_extra/ipsc_gsea/")

g_list3 = c('RBPJ', 'TSC22D1', 'NFIA', 'NR2F1', 'ZFP36L1')

ipsc_markers[ipsc_markers$cluster == 5, ]

FeaturePlot(ipsc, features = g_list1)
FeaturePlot(ipsc, features = g_list2)
FeaturePlot(ipsc, features = g_list3)
ggsave("plot_extra/ipsc/gene_list3.pdf")

OPCs_tdTom1 = readRDS("OPC_ver1_labeled.rds")
OPCs_tdTom1$source = "tdTom"
anchors1 = FindIntegrationAnchors(
        object.list = list(
                tdTom = OPCs_tdTom1,
                ipsc = ipsc
        ))
integrated1 <- IntegrateData(anchorset = anchors1,
                             weight.reduction = "pca",
                            dims = 1:30)
integrated1 <-
        ScaleData(integrated1,
                  features = VariableFeatures(integrated1),
                  verbose = FALSE)
integrated1 <-
        RunPCA(integrated1, features = VariableFeatures(integrated1), npcs = 100)
ElbowPlot(integrated1, ndims = 30)
integrated1 <- RunUMAP(integrated1, reduction = "pca", dims = 1:15)
library(ggplot2)
g1 = UMAPPlot(integrated1, group.by = "source", pt.size = 1.0) +
        theme(legend.text = element_text(size = 20), legend.position = "top")
g2 = UMAPPlot(integrated1[, integrated1$source == "ipsc"], pt.size = 1.0) +
        theme(legend.text = element_text(size = 20), legend.position = "top")
integrated1$cl
cl_col = readr::read_csv("cl_color.csv")
g3 = UMAPPlot(integrated1[, integrated1$source == "tdTom"], pt.size = 1.0) +
        theme(legend.text = element_text(size = 20), legend.position = "top") + scale_color_manual(values = cl_col$Color)
g1 + g2 + g3
ggplotly(g3)

genes_use = intersect(rownames(ipsc), rownames(integrated))
anchors0 = FindTransferAnchors(
        reference = integrated,
        query = ipsc,
        reduction = "cca",
        features = genes_use,
        k.anchor = 5,
        k.score = 30,
        k.filter = 200,
        dims = 1:30
)
predictions0 <- TransferData(
        anchorset = anchors0,
        refdata = integrated$celltype,
        k.weight = 20,
        weight.reduction = "cca",
        dims = 1:20)
ipsc <- AddMetaData(ipsc,
                    metadata = predictions0)
ggplotly(UMAPPlot(ipsc, group.by = "predicted.id", pt.size = 1.0) + 
        theme(legend.text = element_text(size = 20)))
# g_x + g_y
# ggsave("plot_extra/UMAP_predicted_cluster.pdf")

# coembedding
integrated = RenameAssays(integrated, RNA = 'rna')
integrated = RenameAssays(integrated, integrated = 'RNA')
opc_transfer = TransferData(
        anchorset = anchors0,
        refdata = GetAssayData(integrated, assay = "RNA", slot = "data")[genes_use,],
        k.weight = 20,
        weight.reduction = "cca",
        # weight.reduction = OPCs_tdTom1$pca,
        dims = 1:20
)
ipsc = RenameAssays(ipsc, RNA = "rna1")
ipsc[["RNA"]] = opc_transfer

integrated$source = "Reference"
ipsc$source = "ipsc"
coembed <- merge(integrated, y = ipsc)
DefaultAssay(coembed)
coembed <-
        ScaleData(coembed, features = genes_use, do.scale = FALSE)
coembed <-
        RunPCA(
                coembed,
                assay = "RNA",
                features = genes_use,
                verbose = TRUE,
                npcs = 100
        )
# ElbowPlot(coembed, ndims = 25)
is_na = is.na(coembed$celltype)
coembed$celltype[is_na] = coembed$predicted.id[is_na]
coembed <- RunUMAP(coembed, dims = 1:20, repulsion.strength = 1)
g1 = UMAPPlot(coembed, group.by = "source", pt.size = 3.0) + 
        theme(legend.position = "top", legend.text = element_text(size = 30))
g2 = UMAPPlot(coembed, group.by = "celltype", pt.size = 3.0) + 
        scale_color_manual(values = cl_ident$color) + theme(legend.position = "top") +
        theme(legend.position = "top", legend.text = element_text(size = 25),
              plot.margin = margin(r = 4, unit = "cm"))
g1 + g2

DefaultAssay(coembed) = "rna1"
FeaturePlot(coembed[, coembed$source == "ipsc"], features = g_list2[g_list2 %in% genes_use])
FeaturePlot(coembed[, coembed$source == "ipsc"], features = g_list1[g_list1 %in% genes_use])
DefaultAssay(ipsc) = "rna1"

library(dplyr)
top_genes <- ipsc_markers %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
top_genes = top_genes[top_genes$p_val_adj < 0.05, ]
DoHeatmap(ipsc, features = top_genes$gene) + NoLegend()+ theme(plot.margin = margin(2, 1, 1, 1, "cm"))
DoHeatmap(OPCs_tdTom1, features = top_genes$gene, group.colors = col_vec) + NoLegend()+ theme(plot.margin = margin(2, 1, 1, 1, "cm"))

