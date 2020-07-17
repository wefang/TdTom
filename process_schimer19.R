library(Seurat)
library(plotly)
library(readr)
# data obtained from https://cells.ucsc.edu/?ds=ms
schirmer19_mat = Read10X("data/schirmer19/")
schirmer19 = CreateSeuratObject(counts = schirmer19_mat,
                                project = "schirmer19")
meta_df <- read_tsv("data/schirmer19/meta.txt")
schirmer19@meta.data = cbind(schirmer19@meta.data,
                             meta_df[match(colnames(schirmer19), meta_df$cell), -1])
schirmer19 <- subset(schirmer19, subset = nCount_RNA > 2000 & nFeature_RNA > 1000)
schirmer19 = subset(schirmer19, subset = diagnosis == "Control")
# FeatureScatter(schirmer19, "nCount_RNA", "nFeature_RNA")
schirmer19 <- NormalizeData(schirmer19,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
schirmer19 <- FindVariableFeatures(schirmer19,
                                   selection.method = "vst",
                                   nfeatures = 2000)
# write_tsv(schirmer19@meta.data[, "cell_type", drop = F], path = "scvi_dat/schirmer19/labels.tsv")
schirmer19 <- ScaleData(schirmer19)
schirmer19 <- RunPCA(schirmer19, features = VariableFeatures(schirmer19))
ElbowPlot(schirmer19, ndims = 50)
schirmer19 <- RunUMAP(schirmer19, dims = 1:20)
saveRDS(schirmer19, "processed/schirmer19.rds")
UMAPPlot(schirmer19, group.by = "cell_type")

# schirmer19_subset = schirmer19[, schirmer19$cell_type %in% c("Astrocytes",
#                                                              "Endo cells",
#                                                              "OL-Cntl",
#                                                              "OPC")]
# schirmer19_subset <- NormalizeData(schirmer19_subset,
#                                    normalization.method = "LogNormalize",
#                                    scale.factor = 10000)
# schirmer19_subset <- FindVariableFeatures(schirmer19_subset,
#                                           selection.method = "vst",
#                                           nfeatures = 2000)
# schirmer19_subset <- ScaleData(schirmer19_subset)
# schirmer19_subset <- RunPCA(schirmer19_subset, features = VariableFeatures(schirmer19))
# schirmer19_subset <- RunUMAP(schirmer19_subset, dims = 1:10)
# VlnPlot(schirmer19, features = c("BCL11B", "CALB1", "CNP"), group.by = "cell_type")

# ref_var_genes = VariableFeatures(schirmer19_subset)
# opc_var_genes = VariableFeatures(OPCs_tdTom1)
# genes_use = intersect(rownames(OPCs_tdTom1), rownames(schirmer19_subset))
# genes_use = genes_use[genes_use %in% c(ref_var_genes, opc_var_genes)]
# 
# anchors1 = FindTransferAnchors(reference = schirmer19_subset,
#                                query = OPCs_tdTom1,
#                                reduction = "cca",
#                                features = genes_use,
#                                k.anchor = 5,
#                                dims = 1:30)
# # saveRDS(anchors, 'anchors_cca.rds')
# predictions1 <- TransferData(anchorset = anchors1,
#                              refdata = schirmer19_subset$cell_type,
#                              weight.reduction = "cca",
#                              dims = 1:30)
# OPCs_tdTom1 <- AddMetaData(OPCs_tdTom1, metadata = predictions1)
# p2 = ggplotly(UMAPPlot(OPCs_tdTom1, group.by = "predicted.id"))
# p2
# htmlwidgets::saveWidget	(p2, file.path(getwd(), "explore", "schirmer19_labels.html"))
