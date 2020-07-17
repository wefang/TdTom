library(Seurat)
library(plotly)
library(dplyr)
library(ComplexHeatmap)

# normalize_row <- function(x, m = 1) {
#         total = rowSums(x)
#         total[total==0] = 1
#         x / total * m
# }
make_rnk <- function(y, path) {
        # using signed log pval here
        if (!dir.exists(path)) {
                dir.create(path)
        }
        for (i in unique(y$cluster)) {
                x = y[y$cluster == i, ]
                x = x %>% mutate(log_p_signed = log10(p_val_adj) * ifelse(avg_logFC < 0, 1., -1.))
                x = x %>% mutate(stat = -log10(p_val_adj+1e-50) * avg_logFC)
                write_tsv(x = x[order(x$log_p_signed), c("gene", "stat")],
                          path = paste0(path, "/cluster_", i, ".rnk"), col_names = F)
        }
}
OPCs_tdTom1 = readRDS("OPC_ver1_labeled.rds")
OPCs_tdTom1 = FindVariableFeatures(OPCs_tdTom1, nfeatures = 3000)

# cell cycle scoring
g1s_genes = readLines("g1s_genes.txt")
g2m_genes = readLines("g2m_genes.txt")
OPCs_tdTom1 = CellCycleScoring(OPCs_tdTom1, s.features = g1s_genes, g2m.features = g2m_genes)

# FeaturePlot(OPCs_tdTom1, "S.Score", pt.size = 3.0)  + theme(plot.title = element_text(size = 60))
# ggsave(paste0("plots_clean/part2/UMAP_S_score.pdf"))
# FeaturePlot(OPCs_tdTom1, "G2M.Score", pt.size = 3.0) + theme(plot.title = element_text(size = 60))
# ggsave(paste0("plots_clean/part2/UMAP_G2M_score.pdf"))
# UMAPPlot(OPCs_tdTom1, group.by = "Phase", pt.size = 3.0) + ggtitle("Phase") + theme(plot.title = element_text(size = 60))
# ggsave(paste0("plots_clean/part2/UMAP_phase.pdf"))
# 
# FeaturePlot(OPCs_tdTom1, c("S.Score", "G2M.Score"))

# Analysis 1
# regress out cycle score
# OPCs_tdTom1 = ScaleData(OPCs_tdTom1, vars.to.regress =  c("G2M.Score", "S.Score"), features = VariableFeatures(OPCs_tdTom1))
# OPCs_tdTom1 <- RunPCA(OPCs_tdTom1, features = VariableFeatures(OPCs_tdTom1))
# ElbowPlot(OPCs_tdTom1, ndims = 50)
# OPCs_tdTom1 <- RunUMAP(OPCs_tdTom1, dims = 1:30, verbose = FALSE)
# OPCs_tdTom1 <- FindNeighbors(OPCs_tdTom1, dims = 1:30, verbose = FALSE)
# OPCs_tdTom1 <- FindClusters(OPCs_tdTom1, verbose = FALSE)
# cl_map  =as.matrix(table(OPCs_tdTom1$res.1.2, OPCs_tdTom1$seurat_clusters))
# class(cl_map) = "matrix"
# 
# Heatmap(normalize_row(cl_map[as.character(old_cluster_order), ]), cluster_rows = F)
# p = ggplotly(UMAPPlot(OPCs_tdTom1))
# htmlwidgets::saveWidget(p, file.path(getwd(), "cell_cycle/regress/new_clusters.html"))

# OPCs_tdTom1_OPC = OPCs_tdTom1[, !OPCs_tdTom1$seurat_clusters %in% c(2, 3, 8)]
# markers_new = FindAllMarkers(OPCs_tdTom1_OPC, slot = "scale.data", min.pct = 0, logfc.threshold = 0.)
# markers_new_all = FindAllMarkers(OPCs_tdTom1, slot = "scale.data", min.pct = 0, logfc.threshold = 0.)
# 
# top_genes <- markers_new %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
# top_genes = top_genes[top_genes$p_val_adj < 0.05, ]
# top_genes_all <- markers_new_all %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
# top_genes_all = top_genes_all[top_genes_all$p_val_adj < 0.05, ]
# DoHeatmap(OPCs_tdTom1, features = top_genes$gene) + NoLegend()
# DoHeatmap(OPCs_tdTom1, features = top_genes_all$gene) + NoLegend()
# make_rnk(markers_new, path = "cell_cycle/regress/gsea_opc/")
# make_rnk(markers_new_all, path = "cell_cycle/regress/gsea_all/")

# Analysis 2
OPCs_tdTom1$CC.Difference <- OPCs_tdTom1$S.Score - OPCs_tdTom1$G2M.Score
OPCs_tdTom1 <- ScaleData(OPCs_tdTom1, vars.to.regress = "CC.Difference", features = rownames(OPCs_tdTom1))
# Run Analysis on regressed
OPCs_tdTom1 <- RunPCA(OPCs_tdTom1, features = VariableFeatures(OPCs_tdTom1))
# ElbowPlot(OPCs_tdTom1, ndims = 50)
OPCs_tdTom1 <- RunUMAP(OPCs_tdTom1, dims = 1:30, verbose = FALSE)
OPCs_tdTom1 <- FindNeighbors(OPCs_tdTom1, dims = 1:30, verbose = FALSE)
OPCs_tdTom1 <- FindClusters(OPCs_tdTom1, verbose = FALSE)

new_cl_names = c("OPC1",
                 "CyP",
                 "OL",
                 "AS",
                 "OPC2",
                 "OPC3",
                 "OPC4",
                 "Peri",
                 "OPC5")
OPCs_tdTom1$new_cl = factor(new_cl_names[as.numeric(OPCs_tdTom1$seurat_clusters)], levels = c("CyP", paste0("OPC", 1:5), "OL", "AS", "Peri"))
cl_map  =as.matrix(table(OPCs_tdTom1$cl, OPCs_tdTom1$new_cl))
cl_map = cl_map / rowSums(cl_map)
class(cl_map) = "matrix"
Heatmap(cl_map,
        cluster_rows = F, cluster_columns = F,
        column_names_gp = gpar(fontsize = 15),
        row_names_gp = gpar(fontsize = 15),
        name = "Fraction")
UMAPPlot(OPCs_tdTom1, label = T, label.size = 15, pt.size = 3.) + NoLegend()
# htmlwidgets::saveWidget(p, file.path(getwd(), "cell_cycle/regress_diff/new_clusters.html"))

Idents(OPCs_tdTom1) = OPCs_tdTom1$new_cl
OPCs_tdTom1_OPC = OPCs_tdTom1[, !OPCs_tdTom1$seurat_clusters %in% c(2, 3, 7)]

markers_new = FindAllMarkers(OPCs_tdTom1_OPC)
markers_new_all = FindAllMarkers(OPCs_tdTom1)

top_genes <- markers_new %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
top_genes = top_genes[top_genes$p_val_adj < 0.05, ]
top_genes_all <- markers_new_all %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
top_genes_all = top_genes_all[top_genes_all$p_val_adj < 0.05, ]

cl_col = read_csv("cl_color.csv")
col_vec  = cl_col$Color[c(5, 8:12, 7, 2, 13)]
DoHeatmap(OPCs_tdTom1, features = top_genes$gene, group.colors = col_vec, size = 10) + NoLegend() + theme(plot.margin = margin(2, 1, 1, 1, "cm"))
DoHeatmap(OPCs_tdTom1, features = top_genes_all$gene, group.colors = col_vec, size = 10) + NoLegend() + theme(plot.margin = margin(2, 1, 1, 1, "cm"))
UMAPPlot(OPCs_tdTom1, group.by = "new_cl", label = T, label.size = 10, pt.size = 3.) + NoLegend() + scale_color_manual(values = col_vec)
# ggsave("plots_clean/regress_umap.pdf")
# write_csv(markers_new, "regress_markers_opc.csv")
# write_csv(markers_new_all, "regress_markers_all.csv")
# make_rnk(markers_new, path = "cell_cycle/regress_diff/gsea_opc1/")
# make_rnk(markers_new_all, path = "cell_cycle/regress_diff/gsea_all1/")


