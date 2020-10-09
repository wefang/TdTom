library(Seurat)
hESC_tenx = CreateSeuratObject(Read10X("extra_data/hESC_10x/"), min.cells = 10, min.features = 100)
mito.genes <- grep(pattern = "^MT-", x = rownames(hESC_tenx), value = TRUE)
percent.mito <- Matrix::colSums(hESC_tenx@assays$RNA@counts[mito.genes, ])/Matrix::colSums(hESC_tenx@assays$RNA@counts)
hESC_tenx <- AddMetaData(hESC_tenx, metadata = percent.mito, col.name = "percent.mito")
hESC_tenx = hESC_tenx[, hESC_tenx$percent.mito < 0.18]
VlnPlot(hESC_tenx, c("nCount_RNA", "nFeature_RNA", "percent.mito"))

hESC_tenx = FindVariableFeatures(hESC_tenx, nfeatures = 2000)
hESC_tenx = NormalizeData(hESC_tenx)
hESC_tenx = ScaleData(hESC_tenx)
hESC_tenx = RunPCA(hESC_tenx)
hESC_tenx = RunUMAP(hESC_tenx, dims = 1:10)
hESC_tenx = FindNeighbors(hESC_tenx)
hESC_tenx = FindClusters(hESC_tenx)
g_x = UMAPPlot(hESC_tenx, label = T, label.size = 10, pt.size = 3.0, repel = T) + NoLegend() + scale_color_manual(values = cl_col)
g_x

g1s_genes = readLines("g1s_genes.txt")
g2m_genes = readLines("g2m_genes.txt")
hESC_tenx = CellCycleScoring(hESC_tenx, s.features = g1s_genes, g2m.features = g2m_genes)

FeaturePlot(hESC_tenx, features = c("S.Score", "G2M.Score"))

integrated = readRDS("processed/integrated.rds")
genes_use = intersect(rownames(hESC_tenx), rownames(integrated))
anchors0 = FindTransferAnchors(
        reference = integrated,
        query = hESC_tenx,
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
        dims = 1:20
)
hESC_tenx <- AddMetaData(hESC_tenx,
                         metadata = predictions0)
g_y = UMAPPlot(hESC_tenx, group.by = "predicted.id", pt.size = 1.0) +
        scale_color_manual(values = cl_ident$color) +
        theme(legend.text = element_text(size = 20))
g_x + g_y
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
hESC_tenx = RenameAssays(hESC_tenx, RNA = "rna1")
hESC_tenx[["RNA"]] = opc_transfer

integrated$source = "Reference"
hESC_tenx$source = "D85 tdTom+"
coembed <- merge(integrated, y = hESC_tenx)
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
ElbowPlot(coembed, ndims = 25)
is_na = is.na(coembed$celltype)
coembed$celltype[is_na] = coembed$predicted.id[is_na]
coembed <- RunUMAP(coembed, dims = 1:20, repulsion.strength = 1)
g1 = UMAPPlot(coembed, group.by = "source", pt.size = 3.0) +
        theme(legend.position = "top", legend.text = element_text(size = 30))
g2 = UMAPPlot(coembed, group.by = "celltype", pt.size = 3.0) +
        scale_color_manual(values = cl_ident$color) + theme(legend.position = "top") + theme(legend.position = "top", legend.text = element_text(size = 25)) + theme(plot.margin = margin(r = 4, unit = "cm"))
g1 + g2

p_mat = as.matrix(predictions0[, 2:(ncol(predictions0) - 1)])
p_mat = p_mat[, colMaxs(p_mat) > 0.5]
p_mat_norm = p_mat / rowSums(p_mat)
Heatmap(p_mat_norm, show_row_names = F, name = "Probability",
        left_annotation = rowAnnotation(Line = OPCs_tdTom1$res.1.2),
        row_order = order( OPCs_tdTom1$res.1.2))

cluster_ident = do.call(rbind, lapply(split(data.frame(p_mat), hESC_tenx$seurat_clusters), colMeans))
colnames(cluster_ident) = stringr::str_match(colnames(cluster_ident), "prediction\\.score\\.(.*)")[, 2]
rownames(cluster_ident) = cl_labels[match(rownames(cluster_ident), cl_order)]
Heatmap(
        cluster_ident,
        name = "Prob",
        cluster_rows = T,
        column_names_gp = gpar(fontsize = 12),
        row_names_gp = gpar(fontsize = 12),
        right_annotation = rowAnnotation(
                Ident = rownames(cluster_ident),
                col = list(Ident = cl_col),
                show_legend  = F
        )
)
cl_order = c(12, 9, 13, 2, 11, 4, 1, 0, 3, 5, 10, 7, 14, 6, 8)
cl_labels = c("Peri1", "Peri2", "Peri3", "CyP1", "CyP2", "OPC1", "OPC2", "OPC3", "OPC4", "OL1", "OL2", "AS1", "AS2", "AS3", "AS4")
cl_col = c(brewer.pal(6, "Greens")[4:6],
           brewer.pal(4, "Purples")[3:4],
           brewer.pal(6, "RdPu")[2:5],
           brewer.pal(4, "Blues")[3:4],
           brewer.pal(6, "Oranges")[3:6])
names(cl_col) = cl_labels
Idents(hESC_tenx) <- factor(x = hESC_tenx$seurat_clusters, levels = cl_order, labels = cl_labels)

markers = FindAllMarkers(hESC_tenx)
readr::write_csv(markers, "plot_extra/hESC_10x/markers_all.csv")
make_rnk <- function(y, path) {
        # using signed log pval here
        if (!dir.exists(path)) {
                dir.create(path)
        }
        for (i in unique(y$cluster)) {
                x = y[y$cluster == i, ]
                x = x %>% mutate(log_p_signed = log10(p_val_adj) * ifelse(avg_logFC < 0, 1., -1.))
                x = x %>% mutate(stat = -log10(p_val_adj+1e-50) * avg_logFC)
                readr::write_tsv(x = x[order(x$log_p_signed), c("gene", "stat")],
                                 path = paste0(path, "/cluster_", i, ".rnk"), col_names = F)
        }
}
make_rnk(markers, path = "plot_extra/hESC_10x/gsea_all/")

# markers below
hESC_tenx_OPC = hESC_tenx[, hESC_tenx$seurat_clusters %in% c(0, 1, 2, 3, 4, 11)]
markers = FindAllMarkers(hESC_tenx_OPC, logfc.threshold = 0.15)
# make_rnk(markers, path = "plot_extra/hESC_10x/gsea_opcs/")
# readr::write_csv(markers, "plot_extra/markers_OPC.csv")

library(dplyr)
top_genes <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = -p_val_adj)
# top_genes <- markers_no_cycle %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
top_genes = top_genes[top_genes$p_val_adj < 0.05, ]

DoHeatmap(hESC_tenx, features = top_genes$gene, size = 7, angle = 60, group.colors = cl_col) + NoLegend()+ theme(plot.margin = margin(2, 1, 1, 1, "cm"))
# ggsave("plot_extra/hESC_10x/markers_OPC.pdf")
ggsave("plot_extra/hESC_10x/markers_all.pdf")


