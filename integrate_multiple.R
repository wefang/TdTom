library(Seurat)
library(plotly)
library(ComplexHeatmap)
library(RColorBrewer)
library(matrixStats)

zhong = readRDS("processed/zhong.rds")
schirmer19 = readRDS("processed/schirmer19.rds")
jakel = readRDS("processed/jakel.rds")
OPCs_tdTom1 = readRDS("OPC_ver1_labeled.rds")
schirmer19_subset = schirmer19[, schirmer19$cell_type %in% c("Astrocytes",
                                                             "Endo cells",
                                                             "OL-Cntl",
                                                             "OPC")]
jakel_subset = jakel[,!jakel$Celltypes %in% c(paste0("Neuron", 1:5), "Microglia_Macrophages")]
zhong_subset = zhong[, zhong$CT %in% c("Stem cell", "Oligodendrocyte progenitor cell", "Astrocyte")]
# integrate reference
anchors <-
  FindIntegrationAnchors(
    object.list = list(
      schirmer19 = schirmer19_subset,
      jakel = jakel_subset,
      zhong = zhong_subset
    ),
    anchor.features = 2000,
    k.filter = 100,
    dims = 1:30
  )
genes_all = intersect(rownames(zhong),
                      intersect(rownames(jakel),
                                rownames(schirmer19)))
integrated <- IntegrateData(anchorset = anchors,
                            dims = 1:30)
VariableFeatures(integrated) = anchors@anchor.features
integrated <-
  ScaleData(integrated,
            features = VariableFeatures(integrated),
            verbose = FALSE)
integrated <-
  RunPCA(integrated, features = VariableFeatures(integrated), npcs = 100)
ElbowPlot(integrated, ndims = 30)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:15)
cell_type_combine = character(ncol(integrated))
cell_type_combine[!is.na(integrated$cell_type)] = paste0("Schirmer: ", integrated$cell_type[!is.na(integrated$cell_type)])
cell_type_combine[!is.na(integrated$Celltypes)] = paste0("Jakel: ", integrated$Celltypes[!is.na(integrated$Celltypes)])
cell_type_combine[!is.na(integrated$CT)] = paste0("Zhong: ", integrated$CT[!is.na(integrated$CT)])
integrated$cell_type_combine = cell_type_combine
# FeaturePlot(integrated, "PDGFRA")
UMAPPlot(integrated, group.by = "cell_type_combine", pt.size = 3.0) +
  theme(legend.text = element_text(size = 20), legend.position = "top")
# htmlwidgets::saveWidget(p_p, paste0(getwd(), "/explore/integrated_umap.html"))

integrated = FindNeighbors(integrated, dims = 1:30)
integrated = FindClusters(integrated)
UMAPPlot(integrated)
table(integrated$cell_type_combine, integrated$seurat_clusters)
cl_ident = data.frame(
  cluster = 0:13,
  celltype = c(
    "OPC",
    "Oligo1",
    "Astrocyte1",
    "Oligo2",
    "Oligo3",
    "Oligo4",
    "Oligo5",
    "Oligo6",
    "Astrocyte2",
    "Neuro Progenitors",
    "Pericyte/Endo",
    "COPs-ImOligo",
    "OPC-COPs",
    "Astrocyte3"
  )
)
col_def = list(
  OPC = brewer.pal(5, "PuRd")[3:5],
  Oligo = brewer.pal(9, "Blues")[9:3],
  Astrocyte = brewer.pal(5, "OrRd")[3:5],
  Pericyte = brewer.pal(4, "BuGn")[4],
  NP = brewer.pal(3, "BrBG")[1]
)
cl_ident$color = c(
  col_def[["OPC"]][1],
  col_def[["Oligo"]][1],
  col_def[["Astrocyte"]][1],
  col_def[["Oligo"]][c(3:6, 2)],
  col_def[["Astrocyte"]][2],
  col_def[["NP"]],
  col_def[["Pericyte"]],
  col_def[["Oligo"]][7],
  col_def[["OPC"]][2],
  col_def[["Astrocyte"]][3]
)
cl_ident = cl_ident[order(cl_ident$celltype),]
cl_ident$line = c(rep("AL", 3),
                  "OL",
                  "NP",
                  rep("OL", 6),
                  rep("OPC", 2),
                  "Peri")
integrated$celltype = cl_ident[match(integrated$seurat_clusters, cl_ident$cluster), "celltype"]
cl_ident$line = c(rep("AL", 3),
                  "OL",
                  "NP",
                  rep("OL", 6),
                  rep("OPC", 2),
                  "Peri")
cl_ident = cl_ident[order(cl_ident$celltype),]
integrated$line = cl_ident[match(integrated$seurat_clusters, cl_ident$cluster), "line"]
UMAPPlot(
  integrated,
  group.by = "celltype",
  label = T,
  label.size = 12.5,
  repel  = T,
  pt.size = 2.0
) + scale_color_manual(values = cl_ident$color) + NoLegend()
# ggsave("plots_clean/F4.pdf")

ref_var_genes = VariableFeatures(integrated)
opc_var_genes = VariableFeatures(OPCs_tdTom1)
genes_use = intersect(rownames(OPCs_tdTom1), rownames(integrated))
length(genes_use)
genes_use = genes_use[genes_use %in% c(ref_var_genes, opc_var_genes)]
anchors0 = FindTransferAnchors(
  reference = integrated,
  query = OPCs_tdTom1,
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
OPCs_tdTom1 <- AddMetaData(OPCs_tdTom1,
                           metadata = predictions0)
UMAPPlot(OPCs_tdTom1, group.by = "predicted.id", pt.size = 3.0) + 
  scale_color_manual(values = cl_ident$color) + 
  theme(legend.text = element_text(size = 20))
# ggsave("plots_clean/UMAP_predicted_label.pdf")

# coembedding
opc_transfer = TransferData(
  anchorset = anchors0,
  refdata = GetAssayData(integrated, assay = "RNA", slot = "data")[genes_use,],
  k.weight = 20,
  weight.reduction = "cca",
  # weight.reduction = OPCs_tdTom1$pca,
  dims = 1:20
)
OPCs_tdTom1[["RNA"]] = opc_transfer
integrated = RenameAssays(integrated, RNA = 'rna')
integrated = RenameAssays(integrated, integrated = 'RNA')

integrated$source = "Reference"
OPCs_tdTom1$source = "tdTom"
coembed <- merge(integrated, y = OPCs_tdTom1)
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
  scale_color_manual(values = cl_ident$color) + theme(legend.position = "top") + theme(legend.position = "top", legend.text = element_text(size = 25)) + theme(plot.margin = margin(r = 4, unit = "cm"))
g1 + g2
# ggsave("plots_clean//F5.pdf")

cl_label = read_excel("cluster_labels.xlsx")
cl_col = read_csv("cl_color.csv")
cl_col_vec = cl_col$Color
names(cl_col_vec) = cl_col$Ident

# assign transferred probability
p_mat = as.matrix(predictions0[, 2:(ncol(predictions0) - 1)])
p_mat = p_mat[, colMaxs(p_mat) > 0.5]
# p_mat_norm = p_mat / rowSums(p_mat)
# Heatmap(p_mat_norm, show_row_names = F, name = "Probability",
#         left_annotation = rowAnnotation(Line = OPCs_tdTom1$res.1.2),
#         row_order = order( OPCs_tdTom1$res.1.2))

cluster_ident = do.call(rbind, lapply(split(data.frame(p_mat), OPCs_tdTom1$cl), colMeans))
colnames(cluster_ident) = stringr::str_match(colnames(cluster_ident), "prediction\\.score\\.(.*)")[, 2]
Heatmap(
  cluster_ident[as.character(cl_label[2,]),],
  name = "Prob",
  cluster_rows = F,
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  right_annotation = rowAnnotation(
    Ident = as.character(cl_label[2,]),
    col = list(Ident = cl_col_vec),
    show_legend  = F
  )
)

# assign transferred line category probability
predictions0 <- TransferData(
  anchorset = anchors0,
  refdata = integrated$line,
  k.weight = 20,
  weight.reduction = "cca",
  dims = 1:20
)
OPCs_tdTom1 <- AddMetaData(OPCs_tdTom1,
                           metadata = predictions0)
p_mat = as.matrix(predictions0[, 2:(ncol(predictions0) - 1)])
p_mat = p_mat[, colMaxs(p_mat) > 0.5]
cluster_ident = do.call(rbind, lapply(split(data.frame(p_mat), OPCs_tdTom1$cl), colMeans))
colnames(cluster_ident) = stringr::str_match(colnames(cluster_ident), "prediction\\.score\\.(.*)")[, 2]
Heatmap(
  cluster_ident[as.character(cl_label[2,]),],
  name = "Prob",
  cluster_rows = F,
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  right_annotation = rowAnnotation(
    Ident = as.character(cl_label[2,]),
    col = list(Ident = cl_col_vec),
    show_legend  = F
  )
)
