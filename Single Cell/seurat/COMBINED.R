library(dplyr)
library(Seurat)
library(patchwork)

pw41.data <- read.table("~/combined/GSM4483741_PW029-701.cts.txt.gz", header = TRUE, sep = "\t")
pw42.data <- read.table("~/combined/GSM4483742_PW029-702.cts.txt.gz", header = TRUE, sep = "\t")
pw43.data <- read.table("~/combined/GSM4483743_PW030-701.cts.txt.gz", header = TRUE, sep = "\t")
pw44.data <- read.table("~/combined/GSM4483744_PW030-702.cts.txt.gz", header = TRUE, sep = "\t")
pw45.data <- read.table("~/combined/GSM4483745_PW030-703.cts.txt.gz", header = TRUE, sep = "\t")
pw46.data <- read.table("~/combined/GSM4483746_PW030-704.cts.txt.gz", header = TRUE, sep = "\t")
pw47.data <- read.table("~/combined/GSM4483747_PW030-705.cts.txt.gz", header = TRUE, sep = "\t")
pw48.data <- read.table("~/combined/GSM4483748_PW030-706.cts.txt.gz", header = TRUE, sep = "\t")
pw49.data <- read.table("~/combined/GSM4483749_PW030-707.cts.txt.gz", header = TRUE, sep = "\t")
pw50.data <- read.table("~/combined/GSM4483750_PW030-708.cts.txt.gz", header = TRUE, sep = "\t")
pw51.data <- read.table("~/combined/GSM4483751_PW031-712.cts.txt.gz", header = TRUE, sep = "\t")
pw52.data <- read.table("~/combined/GSM4483752_PW032-701.cts.txt.gz", header = TRUE, sep = "\t")
pw53.data <- read.table("~/combined/GSM4483753_PW032-702.cts.txt.gz", header = TRUE, sep = "\t")
pw54.data <- read.table("~/combined/GSM4483754_PW032-703.cts.txt.gz", header = TRUE, sep = "\t")
pw55.data <- read.table("~/combined/GSM4483755_PW032-705.cts.txt.gz", header = TRUE, sep = "\t")
pw56.data <- read.table("~/combined/GSM4483756_PW032-709.cts.txt.gz", header = TRUE, sep = "\t")
pw57.data <- read.table("~/combined/GSM4483757_PW032-711.cts.txt.gz", header = TRUE, sep = "\t")
pw58.data <- read.table("~/combined/GSM4483758_PW032-712.cts.txt.gz", header = TRUE, sep = "\t")
pw59.data <- read.table("~/combined/GSM4483759_PW034-701.cts.txt.gz", header = TRUE, sep = "\t")
pw60.data <- read.table("~/combined/GSM4483760_PW034-702.cts.txt.gz", header = TRUE, sep = "\t")
pw61.data <- read.table("~/combined/GSM4483761_PW034-703.cts.txt.gz", header = TRUE, sep = "\t")
pw62.data <- read.table("~/combined/GSM4483762_PW034-705.cts.txt.gz", header = TRUE, sep = "\t")
pw63.data <- read.table("~/combined/GSM4483763_PW036-701.cts.txt.gz", header = TRUE, sep = "\t")
pw64.data <- read.table("~/combined/GSM4483764_PW036-702.cts.txt.gz", header = TRUE, sep = "\t")
pw65.data <- read.table("~/combined/GSM4483765_PW036-703.cts.txt.gz", header = TRUE, sep = "\t")
pw66.data <- read.table("~/combined/GSM4483766_PW036-705.cts.txt.gz", header = TRUE, sep = "\t")
pw67.data <- read.table("~/combined/GSM4483767_PW040-701.cts.txt.gz", header = TRUE, sep = "\t")
pw68.data <- read.table("~/combined/GSM4483768_PW040-703.cts.txt.gz", header = TRUE, sep = "\t")
pw69.data <- read.table("~/combined/GSM4483769_PW040-705.cts.txt.gz", header = TRUE, sep = "\t")   
pw70.data <- read.table("~/combined/GSM4483770_PW040-707.cts.txt.gz", header = TRUE, sep = "\t")
pw71.data <- read.table("~/combined/GSM4483771_PW040-709.cts.txt.gz", header = TRUE, sep = "\t")
pw72.data <- read.table("~/combined/GSM4483772_PW040-711.cts.txt.gz", header = TRUE, sep = "\t")
pw73.data <- read.table("~/combined/GSM4483773_PW040-712.cts.txt.gz", header = TRUE, sep = "\t")
pw418.data <- read.table("~/combined/GSM5073418_PW051702.cts.txt.gz", header = TRUE, sep = "\t")
pw419.data <- read.table("~/combined/GSM5073419_PW051704.cts.txt.gz", header = TRUE, sep = "\t")
pw420.data <- read.table("~/combined/GSM5073420_PW051708.cts.txt.gz", header = TRUE, sep = "\t")
pw421.data <- read.table("~/combined/GSM5073421_PW052703.cts.txt.gz", header = TRUE, sep = "\t")
pw422.data <- read.table("~/combined/GSM5073422_PW052705.cts.txt.gz", header = TRUE, sep = "\t")
pw423.data <- read.table("~/combined/GSM5073423_PW052706.cts.txt.gz", header = TRUE, sep = "\t")
pw424.data <- read.table("~/combined/GSM5073424_PW052709.cts.txt.gz", header = TRUE, sep = "\t")
pw425.data <- read.table("~/combined/GSM5073425_PW053707.cts.txt.gz", header = TRUE, sep = "\t")
pw426.data <- read.table("~/combined/GSM5073426_PW053710.cts.txt.gz", header = TRUE, sep = "\t")
pw427.data <- read.table("~/combined/GSM5073427_PW053711.cts.txt.gz", header = TRUE, sep = "\t")


pw41 <- CreateSeuratObject(counts = pw41.data, min.cells = 100, min.features = 500)
pw42 <- CreateSeuratObject(counts = pw42.data, min.cells = 100, min.features = 500)
pw43 <- CreateSeuratObject(counts = pw43.data, min.cells = 100, min.features = 500)
pw44 <- CreateSeuratObject(counts = pw44.data, min.cells = 100, min.features = 500)
pw45 <- CreateSeuratObject(counts = pw45.data, min.cells = 100, min.features = 500)
pw46 <- CreateSeuratObject(counts = pw46.data, min.cells = 100, min.features = 500)
 pw47 <- CreateSeuratObject(counts = pw47.data, min.cells = 100, min.features = 500)
pw48 <- CreateSeuratObject(counts = pw48.data, min.cells = 100, min.features = 500)
pw49 <- CreateSeuratObject(counts = pw49.data, min.cells = 100, min.features = 500)
pw50 <- CreateSeuratObject(counts = pw50.data, min.cells = 100, min.features = 500)
pw51 <- CreateSeuratObject(counts = pw51.data, min.cells = 100, min.features = 500)
pw52 <- CreateSeuratObject(counts = pw52.data, min.cells = 100, min.features = 500)
pw53 <- CreateSeuratObject(counts = pw53.data, min.cells = 100, min.features = 500)
pw54 <- CreateSeuratObject(counts = pw54.data, min.cells = 100, min.features = 500)
pw55 <- CreateSeuratObject(counts = pw55.data, min.cells = 100, min.features = 500)
pw56 <- CreateSeuratObject(counts = pw56.data, min.cells = 100, min.features = 500)
pw57 <- CreateSeuratObject(counts = pw57.data, min.cells = 100, min.features = 500)
pw58 <- CreateSeuratObject(counts = pw58.data, min.cells = 100, min.features = 500)
pw59 <- CreateSeuratObject(counts = pw59.data, min.cells = 100, min.features = 500)
pw60 <- CreateSeuratObject(counts = pw60.data, min.cells = 100, min.features = 500)
pw61 <- CreateSeuratObject(counts = pw61.data, min.cells = 100, min.features = 500)
pw62 <- CreateSeuratObject(counts = pw62.data, min.cells = 100, min.features = 500)
pw63 <- CreateSeuratObject(counts = pw63.data, min.cells = 100, min.features = 500)
pw64 <- CreateSeuratObject(counts = pw64.data, min.cells = 100, min.features = 500)
pw65 <- CreateSeuratObject(counts = pw65.data, min.cells = 100, min.features = 500)
pw66 <- CreateSeuratObject(counts = pw66.data, min.cells = 100, min.features = 500)
pw67 <- CreateSeuratObject(counts = pw67.data, min.cells = 100, min.features = 500)
pw68 <- CreateSeuratObject(counts = pw68.data, min.cells = 100, min.features = 500)
pw69 <- CreateSeuratObject(counts = pw69.data, min.cells = 100, min.features = 500)
pw70 <- CreateSeuratObject(counts = pw70.data, min.cells = 100, min.features = 500)
pw71 <- CreateSeuratObject(counts = pw71.data, min.cells = 100, min.features = 500)
pw72 <- CreateSeuratObject(counts = pw72.data, min.cells = 100, min.features = 500)
pw73 <- CreateSeuratObject(counts = pw73.data, min.cells = 100, min.features = 500)
pw418 <- CreateSeuratObject(counts = pw418.data, min.cells = 100, min.features = 500)
pw419 <- CreateSeuratObject(counts = pw419.data, min.cells = 100, min.features = 500)
pw420 <- CreateSeuratObject(counts = pw420.data, min.cells = 100, min.features = 500)
pw421 <- CreateSeuratObject(counts = pw421.data, min.cells = 100, min.features = 500)
pw422 <- CreateSeuratObject(counts = pw422.data, min.cells = 100, min.features = 500)
pw423 <- CreateSeuratObject(counts = pw423.data, min.cells = 100, min.features = 500)
pw424 <- CreateSeuratObject(counts = pw424.data, min.cells = 100, min.features = 500)
pw425 <- CreateSeuratObject(counts = pw425.data, min.cells = 100, min.features = 500)
pw426 <- CreateSeuratObject(counts = pw426.data, min.cells = 100, min.features = 500)
pw427 <- CreateSeuratObject(counts = pw427.data, min.cells = 100, min.features = 500)

GSE148842 <- merge(x= pw41, y = c(pw42, pw43, pw44, pw45, pw46, pw47, pw48, pw49, pw50, pw51, pw52, pw53, pw54, pw55, pw56, pw57, pw58, pw59, pw60, pw61, pw62, pw63, pw64, pw65, pw66, pw67, pw68, pw69, pw70, pw71, pw72, pw73, pw418, pw419, pw420, pw421, pw422, pw423, pw424, pw425, pw426, pw427), 
                   add.cell.ids = c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "418", "419", "420", "421", "422", "423", "424", "425", "426", "427"), 
                   project = "GSE148842",
) 

rm(pw41)
rm(pw42)
rm(pw43)
rm(pw44)
rm(pw45)
rm(pw46)
rm(pw47)
rm(pw48)
rm(pw49)
rm(pw50)
rm(pw51)
rm(pw52)
rm(pw53)
rm(pw54)
rm(pw55)
rm(pw56)
rm(pw57)
rm(pw58)
rm(pw59)
rm(pw60)
rm(pw61)
rm(pw62)
rm(pw63)
rm(pw64)
rm(pw65)
rm(pw66)
rm(pw67)
rm(pw68)
rm(pw69)
rm(pw70)
rm(pw71)
rm(pw72)
rm(pw73)
rm(pw418)
rm(pw419)
rm(pw420)
rm(pw421)
rm(pw422)
rm(pw423)
rm(pw424)
rm(pw425)
rm(pw427)
rm(pw426)
# Load thrm(pw420)e PBMC dataset
BRBMET2.data <- Read10X(data.dir = "324 BRBMET2")
BRBMET3.data <- Read10X(data.dir = "325 BRBMET3")
BRBMET87.data <- Read10X(data.dir = "326 BRBMET87")
LUBMET7.data <- Read10X(data.dir = "327 LUBMET7")
LUBMET1.data <- Read10X(data.dir = "328 LUBMET1")

# Initialize the Seurat object with the raw (non-normalized data).
BRBMET2.data <- CreateSeuratObject(counts = BRBMET2.data, project = "GSE234832", min.cells = 3, min.features = 500)
BRBMET3.data <- CreateSeuratObject(counts = BRBMET3.data, project = "GSE234832", min.cells = 3, min.features = 500)
BRBMET87.data <- CreateSeuratObject(counts = BRBMET87.data, project = "GSE234832", min.cells = 3, min.features = 500)
LUBMET7.data <- CreateSeuratObject(counts = LUBMET7.data, project = "GSE234832", min.cells = 3, min.features = 500)
LUBMET1.data <- CreateSeuratObject(counts = LUBMET1.data, project = "GSE234832", min.cells = 3, min.features = 500)

GSE234832.big <- merge(BRBMET2.data, y = c(BRBMET3.data, BRBMET87.data, LUBMET1.data, LUBMET7.data), add.cell.ids = c("BRBMET2", "BRBMET3", "BRBMET87", "LUBMET1", "LUBMET7"), project = "GSE234832")
GSE234832.big

COMBINED <- merge(GSE234832.big, y = c(GSE148842), add.cell.ids = c("GSE234832", "GSE148842"), project = "GBM",
                  merge.data = TRUE)
saveRDS(COMBINED,
        file = 'COMBINED.rds',
        compress = F)
GetAssayData(COMBINED)[1:10, 1:15]

rm(MERGED)
rm(BRBMET2.data)
rm(BRBMET3.data)
rm(BRBMET87.data)
rm(GSE234832.big)
rm(LUBMET1.data)
rm(LUBMET7.data)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
COMBINED[["percent.mt"]] <- PercentageFeatureSet(COMBINED, pattern = "^mt-")

# Show QC metrics for the first 5 cells
head(COMBINED@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(COMBINED, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(COMBINED, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(COMBINED, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

COMBINED <- NormalizeData(COMBINED, normalization.method = "LogNormalize", scale.factor = 10000)
COMBINED <- FindVariableFeatures(COMBINED, selection.method = "vst", nfeatures = 5000)

# Assuming COMBINED is a data frame or data frame-like object
# Create 'top10' as a subset of the data
top10 <- COMBINED[1:10, ]

# Now, you can use 'top10' in the LabelPoints function
plot1 <- VariableFeaturePlot(COMBINED)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
class(top10)

all.genes <- rownames(COMBINED)
COMBINED <- ScaleData(COMBINED)
COMBINED <- RunPCA(COMBINED, features = VariableFeatures(object = COMBINED))

# Examine and visualize PCA results a few different ways
print(COMBINED[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(COMBINED, dims = 1:2, reduction = "pca")

DimPlot(COMBINED, reduction = "pca")
DimHeatmap(COMBINED, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(COMBINED, dims = 1:15, cells = 500, balanced = TRUE)


COMBINED <- FindNeighbors(COMBINED, dims = 1:10)
COMBINED <- FindClusters(COMBINED, resolution = 0.2)
# Look at cluster IDs of the first 5 cells
head(Idents(COMBINED), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
COMBINED <- RunTSNE(COMBINED, dims = 1:10)
COMBINED <- RunUMAP(COMBINED, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(COMBINED, reduction = "umap")
DimPlot(COMBINED, reduction = "umap", label = TRUE)
#to make tsne plot
DimPlot(COMBINED, reduction = "tsne")
DimPlot(COMBINED, reduction = "tsne", label = TRUE)

write.table(COMBINED, file = "combined.csv")

