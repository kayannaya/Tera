library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

metastasis <- readRDS("metastasis.rds")
View(metastasis@meta.data)
DimPlot(metastasis, reduction = 'umap')
DimPlot(metastasis, reduction = 'umap', label= TRUE)
ref.data <- BlueprintEncodeData(ensembl = TRUE)

predictions <- SingleR(test = metastasis, assay.type.test = 1,
                       ref = ref.data, labels = ref.data$label.main)
