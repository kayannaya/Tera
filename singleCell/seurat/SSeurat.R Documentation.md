---
title: SSeurat.R Documentation

---

### Tera's Projects
## Single Cell Analysis

# SSeurat.R Documentation
This script preprocess Gene Expression Omnibus (GEO) counts data and performs a single cell exploratory data analysis.

## Data that you need
Counts file for the desired samples from a GEO project in a `.txt` or `.txt.gz` file format. You can search for the project that has your desired  in this [website](https://www.ncbi.nlm.nih.gov/geo/).

For example, if we want to explore [GSE6891](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6891), scroll down until you see this

![Screenshot 2026-05-21 at 16.26.41](https://hackmd.io/_uploads/BJyCbS3yGx.png)

Then click the "Series Matrix File(s)", which will lead to this page

![Screenshot 2026-05-21 at 16.59.15](https://hackmd.io/_uploads/Bk1Otr2kGg.png)

Then click the "GSE6891_series_matric.txt.gz" to download all samples' files

You can load this data to the script in the next step.

## 1. Data Loading

The data used in this script is obtained from GEO database. Each samples (with .txt format) are loaded to the pipeline through this line

```r
pw41.data <- read.table(
  "/path/to/GSM4483741_PW029-701.cts.txt.gz", 
  header = TRUE, 
  sep = "\t"
)
```

The function `read.table()` reads the .txt file into a data frame, while the parameter `header = TRUE` shows R that the first row contains column names and `sep = "\t"` specifies that the file is tab-delimited.

*.gz files are automatically decompressed by R*

The above line is repeated for each sample files one by one.

> The output will be a matrix where: rows = genes, columns = cells, and values = raw read counts

## 2. Creating Seurat Object

The matrix needs to be converted into a Seurat object, which is compatible with the upcoming pipeline. It turns the matrix into a special R data structure that stores counts, metadata, and all downstream results together.

```r
pw41 <- CreateSeuratObject(
  counts = pw41.data, 
  min.cells = 100, 
  min.features = 500
)
```

The `counts` retrieves the raw count data from our dataset. `min.cells = 100` only keeps the genes that are detected in *at least* 100 cells, so it removes the very lowly expressed genes. `min.features = 500` only keeps the cells that express at least 500 genes, so it removes empty droplets or very low-quality cells.

## 3. Normalizing Each Sample
Different cells are usually sequenced to different depths, e.g. some cells have more reads than the others. 

```r
pw41 <- NormalizeData(pw41)
```

This process corrects it by applying the function `NormalizeData()`, which applies a log normalization. In depth,

1. Divide each gene's count by the total counts in that cell for every cells available
2. Multiply the result by 10,000 (scaling factor) to bring the values to a comparable range
3. Log transformation od `log(value+1)` is applied to each data points

> The normalized values will be stored in the `data` slot of the Seurat object

## 4. Merging All Samples
As this pipeline uses not only use one sample for each variable, we need to merge the data to explore the project as a whole and obtain a comprehensive understanding about the data.

```r
GSE148842 <- merge(
  x = pw41, 
  y = c(pw42, pw43, pw44, ...pw427), 
  add.cell.ids = c("41", "42", "43", ...,"427"), 
  project = "GSE148842"
)
```

`merge()` function combines multiple Seurat objects that we have created from every samples into a single Seurat object. The `x` variable is the first object that we want to merge, and the `y` variable is a list of all remaining objects.

To ensure that every samples have a distinct value from one another, a prefix of the sample name (e.g. a cell barcode `ACGT` from sample 39 becomes `39_ACGT`) is added to each cell. This prevents the formation of duplicates when the merge happens.

> The output would be an R object that is compressed as and `.Rds` file. `saveRDS()` saves the object to a specific path, and use `readRDS()` to load the data if needed.

## 5. Quality Control

To ensure that the cells that are used in this are alive and have a high quality, a quality control based on mitochondrial genes is conducted in a few steps.

```r
GSE148842[["percent.mt"]] <- PercentageFeatureSet(
  GSE148842, 
  pattern = "^mt-"
)
```
`PercentageFeatureSet()` calculates the percentage of counts coming from a set of genes that matching a pattern. The pattern is defined by a prefix for a certain mitochondrial genes (`^mt-` for mouse and `^MT-` for human).

> High mitochondrial percentage (above 5-10%) suggest that the cell membrane has broken down and the cytoplasmic RNA has leaked out, indicating that the cells are dead.

The result is stored as a new column `percent.mt` in the cell metadata.

```r
head(GSE148842@meta.data, 5)
```

The default columns will show `nFeature_RNA` (genes per cell) `nCount_RNA` (total counts per cell). A violin plot then will be created using the following code.

```r
VlnPlot(GSE148842, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

`features` defines which parameter to plot, and `ncol` arranges the three plots side by side.

Other plots that compare two different parameters by using the `plot()` function.

```
plot1 <- FeatureScatter(GSE148842, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(GSE148842, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

>`plot1` shows cells that has a unusually high mitochondrial content.
>`plot2` shows two cells that are captured together (doublets), which can manipulate the entire data as it will have an abnormally high gene and count numbers.

The next step is to filter the abnormalities that the previous plots shows within a certain cut-off points

```r
GSE148842 <- subset(
  GSE148842, 
  subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5
)
```

The function `subset()` filters out on several conditions.
 1. `nFeature_RNA > 500` removes low-quality cells (too few genes)
 2. `nFeature_RNA < 4000` removes likely doublets (too many genes)
 3. `percent.mt < 5` removes dying cells with high mitochondrial content

> The filtered result should meet all of the above requirements to be able to be pass.

## 6. Feature Selection and Scaling

Variable genes that can distinguish different cell types is important to do to speed up the downstream analysis and reduce noises from housekeeping genes.

```r
GSE148842 <- FindVariableFeatures(
  GSE148842, 
  selection.method = "vst", 
  nfeatures = 4000
)
```

The function `FindVariableFeatures()` selects the most informative genes to differentiates the cell types. In this case, **variance stabilizing transformation (VST)** is used to model the mean-variance relationship and to select genes that has a high variance than a certain threshold. This is defined in the `selection.method = "vst"` argument. `nfeatures = 4000` keeps the top 4000 most variable genes, but this value is adjustable to the desired results.

> The selected genes then are listed and plotted (highlighted in red).

To prepare the data for PCA, a normalization is done once more to eliminate scale differences as PCA is sensitive to it.

```r
GSE148842 <- ScaleData(GSE148842, vars.to.regress = "percent.mt")
```

The `ScaleData()` function centers and scales the gene expression values, which 

1. Substracts the mean expression of each gene across all cells
2. Divides by the standard deviation

> Resulting in a mean = 0 and variance = 1 for each gene.

To remove the effect of mitochondrial percentage from the data (linear regression), `vars.to.regress = "percent.mt"` is used so that it doesn't drive clustering.

## 7. Principal Component Analysis (PCA)

To reduce dimentionality of the data, PCA (`RunPCA()`) is used to compress the thousands of gene dimensions into a smaller number of principal components (PCs). Each PC captures a major axis of variation in the data.

```r
GSE148842 <- RunPCA(
  GSE148842, 
  features = VariableFeatures(object = GSE148842)
)
```

`VizDimLoadings(GSE148842, dims = 1:2, reduction = "pca")` shows which genes that contributes most positively and negatively towards PC1 and PC2. 

Then a heatmap is generated to show the correlation between cells vs. top genes for each PC by using the following code

```r
DimHeatmap(GSE148842, dims = 1:15, cells = 500, balanced = TRUE)
```

## 8. Clustering and Visualization

Clusters are made to visualize the similarities between each cells to each other. Similar cells are grouped or plotted near each other, creating a cluster that defines a certain cell type.

```r
GSE148842 <- FindNeighbors(GSE148842, dims = 1:10)
```
`FindNeighbors()` builds a K-nearest neighbor (KNN) graph in PC space which uses only 10 PCS as defined from `dims = 1:10`.

> Each cell is connected to its most similar neighboring cells in this graph.

Then, the KNN graph is analyzed by using Louvain algorithm to detect the clusters/communities in the KNN.

```r
GSE148842 <- FindClusters(GSE148842, resolution = 0.2)
```

The adjustable argument here is the `resolution` or the granularity. For a *stricter and finer* clusters, increase the value. For a *more generalized and broader* clusters, decrease the value.

> Cluster assignments are stored in `@meta.data` and set as the active identity of each cell

Another dimesion reduction process is done using Uniform Manifold Approximation and Projection to visualize the data into a 2D plot. This methods is a non-linear reduction that preserve both local and global structure.

```
GSE148842 <- RunUMAP(GSE148842, dims = 1:10)
DimPlot(GSE148842, reduction = "umap", label = TRUE)
```

> `DimPlot()` creates a 2D plot where similar cells are placed close to each other, with the cells colored by cluster.

For comparison, another dimesion reduction method is used to create a 2D plot from the previous data.

```r
GSE148842 <- RunTSNE(GSE148842, dims = 1:10)
DimPlot(GSE148842, reduction = "tsne", label = TRUE)
```

`RunTSNE()` runs **t-Distributed Stochastic Neighbor Embedding**, which is better at preserving local structure (finer clusters) than UMAP, but lacks the performance to capture global structure unlike UMAP.

> The output would be two different 2D plots which clusters the cells according to their cell type.

