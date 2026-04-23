library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
BRBMET2.data <- Read10X(data.dir = "324 BRBMET2")
BRBMET3.data <- Read10X(data.dir = "325 BRBMET3")
BRBMET87.data <- Read10X(data.dir = "326 BRBMET87")
LUBMET7.data <- Read10X(data.dir = "327 LUBMET7")
LUBMET1.data <- Read10X(data.dir = "328 LUBMET1")
GSM5645888.data <- Read10X(data.dir = "GSM5645888")
GSM5645889.data <- Read10X(data.dir = "GSM5645889")
GSM5645890.data <- Read10X(data.dir = "GSM5645890")
GSM5645891.data <- Read10X(data.dir = "GSM5645891")
GSM5645892.data <- Read10X(data.dir = "GSM5645892")
GSM5645893.data <- Read10X(data.dir = "GSM5645893")
GSM5645894.data <- Read10X(data.dir = "GSM5645894")
GSM5645895.data <- Read10X(data.dir = "GSM5645895")
GSM5645896.data <- Read10X(data.dir = "GSM5645896")
GSM5645897.data <- Read10X(data.dir = "GSM5645897")
GSM5645898.data <- Read10X(data.dir = "GSM5645898")
GSM5645900.data <- Read10X(data.dir = "GSM5645900")
GSM5645902.data <- Read10X(data.dir = "GSM5645902")
GSM5645904.data <- Read10X(data.dir = "GSM5645904")
GSM5645906.data <- Read10X(data.dir = "GSM5645906")
GSM5645908.data <- Read10X(data.dir = "GSM5645908")

# Initialize the Seurat object with the raw (non-normalized data).
BRBMET2.data <- CreateSeuratObject(counts = BRBMET2.data, project = "GSE234832", min.cells = 3, min.features = 500)
BRBMET3.data <- CreateSeuratObject(counts = BRBMET3.data, project = "GSE234832", min.cells = 3, min.features = 500)
BRBMET87.data <- CreateSeuratObject(counts = BRBMET87.data, project = "GSE234832", min.cells = 3, min.features = 500)
LUBMET7.data <- CreateSeuratObject(counts = LUBMET7.data, project = "GSE234832", min.cells = 3, min.features = 500)
LUBMET1.data <- CreateSeuratObject(counts = LUBMET1.data, project = "GSE234832", min.cells = 3, min.features = 500)

GSE234832.big <- merge(BRBMET2.data, y = c(BRBMET3.data, BRBMET87.data, LUBMET1.data, LUBMET7.data), add.cell.ids = c("BRBMET2", "BRBMET3", "BRBMET87", "LUBMET1", "LUBMET7"), project = "GSE234832")
GSE234832.big

GSM5645888.data <- CreateSeuratObject(counts = GSM5645888.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645889.data <- CreateSeuratObject(counts = GSM5645889.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645890.data <- CreateSeuratObject(counts = GSM5645890.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645891.data <- CreateSeuratObject(counts = GSM5645891.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645892.data <- CreateSeuratObject(counts = GSM5645892.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645893.data <- CreateSeuratObject(counts = GSM5645893.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645894.data <- CreateSeuratObject(counts = GSM5645894.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645895.data <- CreateSeuratObject(counts = GSM5645895.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645896.data <- CreateSeuratObject(counts = GSM5645896.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645897.data <- CreateSeuratObject(counts = GSM5645897.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645898.data <- CreateSeuratObject(counts = GSM5645898.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645900.data <- CreateSeuratObject(counts = GSM5645900.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645902.data <- CreateSeuratObject(counts = GSM5645902.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645904.data <- CreateSeuratObject(counts = GSM5645904.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645906.data <- CreateSeuratObject(counts = GSM5645906.data, project = "GSE186344", min.cells = 3, min.features = 500)
GSM5645908.data <- CreateSeuratObject(counts = GSM5645908.data, project = "GSE186344", min.cells = 3, min.features = 500)

GSE186344.big <- merge(GSM5645888.data, y = c(GSM5645889.data, GSM5645890.data, GSM5645891.data, GSM5645892.data, GSM5645893.data, GSM5645894.data, GSM5645895.data, GSM5645896.data, GSM5645897.data, GSM5645898.data, GSM5645900.data, GSM5645902.data, GSM5645904.data, GSM5645906.data, GSM5645908.data), add.cell.ids = c("GSM5645888", "GSM5645889", "GSM5645890", "GSM5645891", "GSM5645892", "GSM5645893", "GSM5645894", "GSM5645895", "GSM5645896", "GSM5645897", "GSM5645898", "GSM5645900", "GSM5645902", "GSM5645904", "GSM5645906", "GSM5645908"), project = "GSE186344")
GSE186344.big

metastasis <- merge(GSE186344.big, y = c(GSE234832.big), add.cell.ids = c("GSE234832", "GSE186344"), project = "glioblastoma")
metastasis

metastasis[c("CD3D", "TCL1A", "MS4A1"), 1:30]
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
metastasis[["percent.mt"]] <- PercentageFeatureSet(metastasis, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(metastasis@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(metastasis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster = FALSE)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(metastasis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(metastasis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
metastasis <- subset(metastasis, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 40)
metastasis <- NormalizeData(metastasis, normalization.method = "LogNormalize", scale.factor = 10000)
metastasis <- NormalizeData(metastasis)
metastasis <- FindVariableFeatures(metastasis, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(metastasis), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(metastasis)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(metastasis)
metastasis <- ScaleData(metastasis, features = all.genes)
metastasis <- ScaleData(metastasis, vars.to.regress = "percent.mt")
metastasis <- RunPCA(metastasis, features = VariableFeatures(object = metastasis))
# Examine and visualize PCA results a few different ways
print(metastasis[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(metastasis, dims = 1:2, reduction = "pca")
DimPlot(metastasis, reduction = "pca")
DimHeatmap(metastasis, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(metastasis, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
metastasis <- JackStraw(metastasis, num.replicate = 100)
metastasis <- ScoreJackStraw(metastasis, dims = 1:20)
JackStrawPlot(metastasis, dims = 1:15)
ElbowPlot(metastasis)
metastasis <- FindNeighbors(metastasis, dims = 1:10)
metastasis <- FindClusters(metastasis, resolution = 0.1)
# Look at cluster IDs of the first 5 cells
head(Idents(metastasis), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
metastasis <- RunUMAP(metastasis, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(metastasis, reduction = "umap")
DimPlot(metastasis, reduction = "umap", label = TRUE)
metastasis <- RunTSNE(metastasis, dims = 1:10)
DimPlot(metastasis, reduction = "tsne")
DimPlot(metastasis, reduction = "tsne", label = TRUE)
saveRDS(metastasis, file = "metastasis.rds")


