library(Seurat)
library(readr)
library(Matrix)

# data collected from "https://db.cngb.org/HCL/gallery.html?tissue=Fetal-Brain_Zhong"
df = read_csv("data/Fetal-Brain_Zhong/extracted/Fetal-Brain_Zhong_dge.txt")
meta = read_csv("data/Fetal-Brain_Zhong/extracted/Fetal-Brain_Zhong_Anno.csv")
mat = as.matrix(df[, -1])
rownames(mat) = df$X1
mat_sp = Matrix(mat, sparse = T)
x_obj = CreateSeuratObject(counts = mat_sp)
all(meta$X1 == colnames(x_obj))
x_obj@meta.data = cbind(x_obj@meta.data, meta[, -1])

FeatureScatter(x_obj, "nCount_RNA", "nFeature_RNA")
x_obj = FindVariableFeatures(x_obj, selection.method = "vst",
                             nfeatures = 2000, verbose = FALSE)
VariableFeaturePlot(x_obj)
x_obj = NormalizeData(x_obj)
x_obj = ScaleData(x_obj)
x_obj = RunPCA(x_obj)
ElbowPlot(x_obj)
x_obj = RunUMAP(x_obj, dims = 1:15)
UMAPPlot(x_obj, group.by = "CT")
# saveRDS(x_obj, file = "processed/zhong.rds")

# x_obj = FindNeighbors(x_obj)
# x_obj = FindClusters(x_obj, resolution = 0.2)

# OPCs_tdTom1 = readRDS("OPC_ver1.rds")
# anchors = FindTransferAnchors(reference = x_obj, query = OPCs_tdTom1,
#                               # reduction = "cca",
#                               features = genes_use, k.anchor = 5,
#                               dims = 1:30)
# predictions <- TransferData(anchorset = anchors, refdata = x_obj$CT,
#                             # weight.reduction = "cca",
#                             dims = 1:30)
# OPCs_tdTom1 <- AddMetaData(OPCs_tdTom1, metadata = predictions)
# p_pred = ggplotly(UMAPPlot(OPCs_tdTom1, group.by = "predicted.id"))
# htmlwidgets::saveWidget(p_pred, "zhong_labels.html")