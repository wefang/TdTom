library(ggplot2)
library(readr)
library(readxl)
library(dplyr)
library(Seurat)

# this function makes rank file for preranked GSEA
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


# rda file produced by seurat 2.0 old analysis
load('./OPCs_042820.RData')
OPCs_tdTom1 = OPCs_tdTom
OPCs_tdTom1 <- NormalizeData(OPCs_tdTom1, verbose = FALSE)
OPCs_tdTom1 <- ScaleData(OPCs_tdTom1, verbose = FALSE)
OPCs_tdTom1 <- FindVariableFeatures(OPCs_tdTom1, selection.method = "vst",
                                    nfeatures = 2000, verbose = FALSE)
OPCs_tdTom1 = ScaleData(OPCs_tdTom1)
OPCs_tdTom1 = RunPCA(OPCs_tdTom1, verbose = TRUE, npcs = 100)
ElbowPlot(OPCs_tdTom1, ndims = 50)
OPCs_tdTom1 = FindNeighbors(OPCs_tdTom1, dims = 1:20)
OPCs_tdTom1 = FindClusters(OPCs_tdTom1)
OPCs_tdTom1 <- RunUMAP(OPCs_tdTom1, dims = 1:20)

cluster10_sub = read_csv("c10_subclustered_ids.csv")
stopifnot(all(OPCs_tdTom1$res.1.2[match(cluster10_sub$X1, colnames(OPCs_tdTom1))] == "10"))
OPCs_tdTom1$res.1.2[match(cluster10_sub$X1, colnames(OPCs_tdTom1))] = cluster10_sub$idents
cl_label = read_excel("cluster_labels.xlsx")
names(cl_label)[13] = "10.2"
OPCs_tdTom1$cl = factor(as.character(cl_label[2, ])[match(OPCs_tdTom1$res.1.2, names(cl_label))], levels = as.character(cl_label[2, ]))
cl_col = read_csv("cl_color.csv")
cl_col
col_vec  = cl_col$Color[match(as.character(cl_label[2, ]), cl_col$Ident)]

UMAPPlot(OPCs_tdTom1, group.by = "cl", label = T, label.size = 10, pt.size = 3.0, repel = T) + NoLegend() + scale_color_manual(values = col_vec)
# ggsave("plots_clean/F1.pdf")
UMAPPlot(OPCs_tdTom1, group.by = "cl", label = T, label.size = 15, pt.size = 3.0, repel = T) + NoLegend() + scale_color_manual(values = col_vec)
# ggsave("plots_clean/F1_larger.pdf")
UMAPPlot(OPCs_tdTom1, group.by = "res.1.2", label = T, label.size = 15, pt.size = 3.0, repel = T) + NoLegend() + scale_color_manual(values =cl_col$Color[match(as.character(cl_label[2, ]), cl_col$Ident)][order(names(cl_label))] )
# ggsave("plots_clean/F1_supp.pdf")

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

g_list1 = c("PMP2", "PRDX1", "TOP2A", "PCNA", "MKI67", "MCM6", "RTN1", "NNAT", "SOX4", "STMN2", "PTPRZ1", "HMP19", "SCG3", "IFI6", "H19", "ISG15", "HS3ST1", "PHLDA1", "TRIO", "SOX6", "FGF12")
g_list2 = c("PDGFRA", "OLIG1", "OLIG2", "NKX2-2", "SOX10", "TCF7L2", "MBP", "CNP", "PLP1", "ZNF488", "CD9", "SOX9", "NFIA", "GFAP", "AQP4", "KCTD12", "MAG", "MOG", "SPOCK1", "BCAS1", "COL1A1", "COL1A2", "COL3A1", "ACTA2", "S100B")
# for (g in c(g_list1, g_list2)) {
#         FeaturePlot(OPCs_tdTom1, g, pt.size = 3.0) + theme(plot.title = element_text(size = 60))
#         ggsave(paste0("plots_clean/part3/", g, ".pdf"))
# }
UMAPPlot(OPCs_tdTom1, group.by = "age", pt.size = 3.0) + theme(plot.title = element_text(size = 60))
# ggsave(paste0("plots_clean/part2/UMAP_age.pdf"))
FeaturePlot(OPCs_tdTom1, "percent.mito", pt.size = 3.0) + theme(plot.title = element_text(size = 60))
# ggsave(paste0("plots_clean/part2/UMAP_mito.pdf"))

cluster_id = do.call(rbind, list(data.frame(old_cluster = c(0, 1, 3, 5, 6, 7, 11), type = "OPC"),
                                 data.frame(old_cluster = 2, type = "ALC"),
                                 data.frame(old_cluster = c(4, 8), type = "OLLC"),
                                 data.frame(old_cluster = 10, type = "Pericyte-Astrocyte"),
                                 data.frame(old_cluster = 9, type = "NeuroProgenitor")))
cluster_id$old_cluster = as.character(cluster_id$old_cluster)
Idents(OPCs_tdTom1) = OPCs_tdTom1$cl
saveRDS(OPCs_tdTom1, file = "OPC_ver1_labeled.rds")

# OPCs_tdTom1 = readRDS("OPC_ver1_labeled.rds")
cluster10_sub = read_csv("c10_subclustered_ids.csv")
stopifnot(all(OPCs_tdTom1$res.1.2[match(cluster10_sub$X1, colnames(OPCs_tdTom1))] == "10"))
OPCs_tdTom1$res.1.2[match(cluster10_sub$X1, colnames(OPCs_tdTom1))] = cluster10_sub$idents
cl_label = read_excel("cluster_labels.xlsx")
names(cl_label)[13] = "10.2"
OPCs_tdTom1$cl = factor(as.character(cl_label[2, ])[match(OPCs_tdTom1$res.1.2, names(cl_label))], levels = as.character(cl_label[2, ]))
cl_col = read_csv("cl_color.csv")
col_vec  = cl_col$Color[match(as.character(cl_label[2, ]), cl_col$Ident)]

OPCs_tdTom1_OPC = OPCs_tdTom1[, OPCs_tdTom1$res.1.2 %in% cluster_id[cluster_id$type %in% c("OPC"), "old_cluster"]]
markers = FindAllMarkers(OPCs_tdTom1_OPC)
# make_rnk(markers, path = "gsea1_opc/")

top_genes <- markers %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
# top_genes <- markers_no_cycle %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
top_genes = top_genes[top_genes$p_val_adj < 0.05, ]
# Idents(OPCs_tdTom1) <- factor(x = OPCs_tdTom1$res.1.2, levels = cluster_order)
# OPCs_tdTom1$res.1.2 = factor(x = OPCs_tdTom1$res.1.2, levels = cluster_order)
# OPCs_tdTom1 <- ScaleData(OPCs_tdTom1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(OPCs_tdTom1))
OPCs_tdTom1 = ScaleData(OPCs_tdTom1, features = rownames(OPCs_tdTom1))
DoHeatmap(OPCs_tdTom1, features = top_genes$gene, group.colors = col_vec, size = 7, angle = 60) + NoLegend()+ theme(plot.margin = margin(2, 1, 1, 1, "cm"))

# OPCs_tdTom1 = ScaleData(OPCs_tdTom1, features = rownames(OPCs_tdTom1))
# OPCs_tdTom1_OPC = ScaleData(OPCs_tdTom1_OPC, features = rownames(OPCs_tdTom1))

OPCs_tdTom1_ncOPC = OPCs_tdTom1[, OPCs_tdTom1$res.1.2 %in% c(0, 1, 3, 11)]
markers = FindAllMarkers(OPCs_tdTom1_ncOPC)
top_genes <- markers %>% group_by(cluster) %>% top_n(n = 40, wt = -p_val_adj)
top_genes = top_genes[top_genes$p_val_adj < 0.05, ]
DoHeatmap(OPCs_tdTom1, features = top_genes$gene, group.colors = col_vec, size = 7, angle = 60) + NoLegend()+ theme(plot.margin = margin(2, 1, 1, 1, "cm"))

old_markers = read_excel("Supplemental table 2.xlsx", col_names = T, skip = 1)
top_genes_old <- old_markers %>% group_by(cluster) %>% top_n(n = 20, wt = -p_val_adj)
DoHeatmap(OPCs_tdTom1, features = top_genes_old$gene, group.colors = col_vec, size = 7, angle = 60) + NoLegend()+ theme(plot.margin = margin(2, 1, 1, 1, "cm"))

markers_all = FindAllMarkers(OPCs_tdTom1)
# make_rnk(markers_all, path = "gsea1_all/")



