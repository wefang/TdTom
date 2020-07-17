library(readr)
library(Seurat)
library(ggplotly)

# data obtained from GSE118257
cell_names= read_tsv("data/jakel19/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz", n_max = 1)
df =read_tsv("data/jakel19/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz", skip = 1, col_names = F)
count_mat = as.matrix(df[, -1])
colnames(count_mat) = colnames(cell_names)
rownames(count_mat) = df$X1
meta_df = read_tsv("data/jakel19/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz")
colnames(meta_df)[1] = "cell"
colnames(meta_df)[2] = "Detected"

jakel = CreateSeuratObject(counts = count_mat)
jakel@meta.data = cbind(jakel@meta.data,
                        meta_df[match(colnames(jakel), meta_df$cell), -1])

# strip_version <- function(x) {
#         sapply(strsplit(x, "\\."), "[[", 1)
# }
# FeatureScatter(jakel, "nCount_RNA", "nFeature_RNA")
jakel <- subset(jakel, subset = nCount_RNA > 1000 & nFeature_RNA > 500)
jakel = subset(jakel, subset = Condition == "Ctrl")
jakel <- NormalizeData(jakel,
                       normalization.method = "LogNormalize",
                       scale.factor = 5000)
jakel <- FindVariableFeatures(jakel,
                              selection.method = "vst",
                              nfeatures = 2000)
jakel <- ScaleData(jakel)
jakel <- RunPCA(jakel, features = VariableFeatures(jakel))
ElbowPlot(jakel, ndims = 50)
jakel <- RunUMAP(jakel, dims = 1:20)

# VariableFeaturePlot(jakel)
# VizDimLoadings(jakel, dims = 1:2, reduction = "pca")
# ElbowPlot(jakel, ndims = 50)
p = ggplotly(UMAPPlot(jakel, group.by = "Celltypes"))
htmlwidgets::saveWidget(widget = p, file = file.path(getwd(), "explore", "jakel.html"))
# saveRDS(jakel, file = "processed/jakel.rds")

# jakel = readRDS("processed/jakel.rds")
# jakel_subset = jakel[, !jakel$Celltypes %in% c(paste0("Neuron", 1:5), "Microglia_Macrophages")]
# jakel_subset <- FindVariableFeatures(jakel_subset,
#                                      selection.method = "vst",
#                                      nfeatures = 2000)
# jakel_subset <- ScaleData(jakel_subset)
# jakel_subset <- RunPCA(jakel_subset, features = VariableFeatures(jakel_subset))
# ElbowPlot(jakel_subset, ndims = 50)
# jakel_subset <- RunUMAP(jakel_subset, dims = 1:20)
# ref_var_genes = VariableFeatures(jakel_subset)
# opc_var_genes = VariableFeatures(OPCs_tdTom1)
# genes_use = intersect(rownames(OPCs_tdTom1), rownames(jakel))
# genes_use = genes_use[genes_use %in% c(ref_var_genes, opc_var_genes)]
# 
# anchors1 = FindTransferAnchors(reference = jakel_subset, query = OPCs_tdTom1,
#                                reduction = "cca",
#                                features = genes_use,
#                                k.anchor = 10,
#                                dims = 1:30)
# predictions1 <- TransferData(anchorset = anchors1, refdata = jakel_subset@meta.data$Celltypes,
#                              weight.reduction = "cca",
#                              k.weight = 50,
#                              dims = 1:30)
# table(predictions1$predicted.id)
# OPCs_tdTom1 <- AddMetaData(OPCs_tdTom1, metadata = predictions1)
# p2 = ggplotly(UMAPPlot(OPCs_tdTom1, group.by = "predicted.id") + 
#                       scale_color_manual(values=col_mapper(sort(unique(predictions1$predicted.id)), jakel_map)))
# p2
# htmlwidgets::saveWidget(p2, file.path(getwd(), "explore", "jakel_labels.html"))
# 
# opc_transfer = TransferData(anchorset = anchors1,
#                             refdata = GetAssayData(jakel_subset, assay = "RNA", slot = "data")[genes_use, ],
#                             weight.reduction = "cca",
#                             # weight.reduction = OPCs_tdTom1@reductions$umap,
#                             dims = 1:30)
# OPCs_tdTom1[["RNA"]] = opc_transfer
# jakel_subset$source = "Jakel"
# OPCs_tdTom1$source = "tdTom"
# coembed <- merge(jakel_subset, OPCs_tdTom1)
# DefaultAssay(coembed)
# coembed <- ScaleData(coembed, features = genes_use, do.scale = FALSE)
# coembed <- RunPCA(coembed, assay = "RNA", features = genes_use, verbose = TRUE, npcs = 100)
# ElbowPlot(coembed, ndims = 25)
# is_na = is.na(coembed$Celltypes)
# coembed$Celltypes[is_na] = coembed$predicted.id[is_na]
# coembed <- RunUMAP(coembed, dims = 1:30, repulsion.strength = 5)
# g1 = UMAPPlot(coembed, group.by = c("source"))
# g2 = UMAPPlot(coembed, group.by = c("Celltypes")) + 
#         scale_color_manual(values = col_mapper(sort(unique(coembed$Celltypes)),
#                                                jakel_map))
# g1 + g2


