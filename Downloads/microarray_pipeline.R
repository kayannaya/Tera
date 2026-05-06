# configurations
# set your directories here
setwd("/Users/kayannaya/Documents/Work/TEEP/Tera/Downloads")
raw_dir <- "/Users/kayannaya/Documents/Work/TEEP/Tera/Downloads/GSE38832_RAW"

# Install the required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c(
  "limma", "affy", "Biobase", "arrayQualityMetrics",
  "GEOquery", "hgu133plus2cdf", "hgu133plus2probe", "hgu133plus2.db"
)

cran_pkgs <- c(
  "ggplot2", "lattice", "ggpubr", "rstatix",
  "purrr", "tidyselect", "dplyr", "NMF", "magrittr"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(limma)
library(affy)
library(Biobase)
library(arrayQualityMetrics)
library(magrittr)
library(ggplot2)
library(GEOquery)
library(lattice)
library(hgu133plus2cdf)
library(hgu133plus2probe)
library(hgu133plus2.db)
library(ggpubr)
library(rstatix)
library(purrr)
library(tidyselect)
library(dplyr)
library(NMF)

# GSE_getdata <- function(GSE_id, cdf_name){
#   wd <- getwd()
#   getGEOSuppFiles(GDE_id)
#   gds <- getGEO(GDE_id)
#   
#   # Prepare metadata
#   pheno_data <- pData(phenoData(gds))
#   targets <- rownames(pheno_data)
#   
#   # Prepare data
#   untars
#   ReadAffy()
# }

# obtain the data
# pick one GEO projects, perform ctrl f and replace all the previous GEO project name with your new one. 
# the file names will automatically be changed as well
options(timeout = 300)
gds <- getGEO("GSE38832") 
pheno_data <- pData(phenoData(gds$GSE38832_series_matrix.txt.gz))

# create dir and download raw CEL files from GEO
dir.create("~/Affy_microArray/GSE38832", recursive = TRUE, showWarnings = FALSE)
getGEOSuppFiles("GSE38832", baseDir = "~/Affy_microArray")

# create the directory
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# untar the downloaded archive
untar("/Users/kayannaya/Affy_microArray/GSE38832/GSE38832_RAW.tar",
      exdir = raw_dir)

# read the CEL files
data_raw <- ReadAffy(celfile.path = raw_dir,
                     cdfname = "hgu133plus2")

# WRONG column for stage
# pheno_data <- subset(pheno_data, select = characteristics_ch1.4)  

# use characteristics_ch1 for ajcc_stage
pheno_data <- subset(pheno_data, select = characteristics_ch1)
colnames(pheno_data)[1] <- "sample_type"
pheno_data$sample_type <- gsub("ajcc_stage: ", "", pheno_data$sample_type)
pheno_data$sample_type <- factor(pheno_data$sample_type)

# check actual levels first
levels(pheno_data$sample_type)

# then assign correctly based on what you see
levels(pheno_data$sample_type)[levels(pheno_data$sample_type) == "1"] <- "I"
levels(pheno_data$sample_type)[levels(pheno_data$sample_type) == "2"] <- "II"
levels(pheno_data$sample_type)[levels(pheno_data$sample_type) == "3"] <- "III"
levels(pheno_data$sample_type)[levels(pheno_data$sample_type) == "4"] <- "IV"

rownames(pheno_data) <- paste(rownames(pheno_data), ".CEL.gz", sep = "") # add .CEL.gz to sampleID

# Load phenotype data (normal/tumor), set rownames to sampleIDs
meta_data <- pData(phenoData(gds$GSE38832_series_matrix.txt.gz))
rownames(meta_data) <- paste0(rownames(meta_data), ".CEL.gz")

# Read all CEL files in path, with phenotype detail
pData(data_raw) # check if phenotype data is added

# ALT
data_rma <- rma(data_raw)
data_rma_norm = exprs(data_rma)
save(data_raw, data_rma, file = "affyData_GSE38832.rda")
# write.exprs(data_rma, "expression_rma_GSE38832.txt") # gets and writes expression estimates
load("affyData_GSE38832.rda")

# normalize data if needed
# # LOG TRANSFORM
# data_rma_norm_log = log2(data_rma_norm  + 1)
# 
# # QUANTILE
# data_rma_norm_qua = normalize.quantiles(data_rma_norm)

probes_sd = apply(data_rma_norm, 1, sd)
thres = 0.2
var_probes = names(probes_sd[probes_sd > thres])
var_matrix = data_rma_norm[var_probes,]

probes = rownames(var_matrix) # gets probes ids
symbols = data.frame(genes = unlist(
  mget(probes, 
       hgu133plus2SYMBOL, 
       ifnotfound = NA))) # gets gene symbol associated to probe
# genes_probes = cbind(probes = rownames(symbols),symbols)
matrix_out = data.frame(cbind(symbols,var_matrix)) # creates new df with added column of gene symbol
matrix_out$genes[is.na(matrix_out$genes)] = rownames(matrix_out)[is.na(matrix_out$genes)]
matrix_out = matrix_out %>%
  group_by(genes) %>%
  summarise_all(mean) # groups genes with the same name and gets exp mean
matrix_out = data.frame(matrix_out)

# swap the order of this
rownames(matrix_out) = matrix_out$genes
matrix_out = matrix_out[,-1]

# and this
genes = c("TEX36-AS1","AMPK","ULK1","ATG13","ULK2", 
          "ATG7", "ATG5","ATG3","ATG12","ATG16L",
          "LAMP2", "p62","mTORC1","Plekhm1")
genes = alias2Symbol(genes)

test = matrix_out[genes,]


## DIMENSIONAL REDUCTION

# t_data is defined later in your script but used here first
t_data = data.frame(t(matrix_out))
hist(as.matrix(matrix_out), ylab="", las=2, main="Raw data")
boxplot(as.matrix(matrix_out), ylab = expression('Log'[2]~'Read counts'), las = 2, main = "Raw data")

# also rank=15000 is the number of features, not the NMF rank (should be 2-10)
res = nmf(t_data, rank = 3, seed = 123)  # rank = number of metagenes, typically 2-10
W = basis(res)
H = coef(res)

## CORR
t_data = data.frame(t(matrix_out))
corr_mat = cor(t_data)
corrplot::corrplot(corr_mat)

t_data %>% cor_test(vars = "X1556133_s_at", 
                    vars2 = c("X1554112_a_at", "X1569827_at", "X202512_s_at"))

# WRONG - hyphens break formula parsing
# test = lm(formula = "TEX36-AS1" ~ genes, data= t_data)

# CORRECT - use backticks for non-standard names
test = lm(`TEX36-AS1` ~ ., data = t_data[, c("TEX36-AS1", genes)])

mem.maxVSize()