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
#   untar
#   ReadAffy()
# }

gds <- getGEO("GSE388332")
pheno_data <- pData(phenoData(gds$GSE38832_series_matrix.txt.gz))
pheno_data <- subset(pheno_data, select = characteristics_ch1.4)
colnames(pheno_data)[1] = "sample_type"
pheno_data$sample_type = factor(pheno_data$sample_type)

levels(pheno_data$sample_type)[1] = "II"
levels(pheno_data$sample_type)[2] = "III"

rownames(pheno_data) <- paste(rownames(pheno_data), ".CEL.gz", sep = "") # add .CEL.gz to sampleID


# Load phenotype data (normal/tumor), set rownames to sampleIDs
meta_data <- read.table("GSE15471_series_matrix_metadata.txt", header = T, row.names = 1)
rownames(meta_data) <- paste(rownames(meta_data), ".CEL.gz", sep = "") # add .CEL.gz to sampleID

# Read all CEL files in path, with phenotype detail
data_raw <- ReadAffy(celfile.path = "~/Affy_microArray/GSE38832_RAW",
                 cdfname = "hgu133plus2")
pData(data_raw) # check if phenotype data is added

# ALT
data_rma <- rma(data_raw, target= "core")
data_rma_norm = exprs(data_rma)
save(data, data_rma, file = "affyData_GSE38832.rda")
# write.exprs(data_rma, "expression_rma_GSE38832.txt") # gets and writes expression estimates
load("affyData_GSE38832.rda")

# LOG TRANSFORM
data_rma_norm_log = log2(data_rma_norm  + 1)

# QUANTILE
data_rma_norm_qua = normalize.quantiles(data_rma_norm)


probes_sd = apply(data_rma_norm_log, 1, sd)
thres = 0.2
var_probes = names(probes_sd[probes_sd > thres])
var_matrix = data_rma_norm_log[var_probes,]

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

genes = c("TEX36-AS1","AMPK","ULK1","ATG13","ULK2", 
          "ATG7", "ATG5","ATG3","ATG12","ATG16L",
          "LAMP2", "p62","mTORC1","Plekhm1")
genes = alias2Symbol(genes)

test = matrix_out[genes,]


rownames(matrix_out) = matrix_out$genes
matrix_out = matrix_out[,-1]


## DIMENSIONAL REDUCTION
hist(matrix_out,ylab="",las=2,main="Raw data")
boxplot(log2(as.matrix(matrix_out)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Raw data")

res = nmf(t_data, rank = 15000, seed = 123)
W = basis(res)
H = coef(res) 

genes = c("TEX36-AS1","AMPK","ULK1","ATG13","ULK2", 
          "ATG7", "ATG5","ATG3","ATG12","ATG16L",
          "LAMP2", "p62","mTORC1","Plekhm1")
genes = alias2Symbol(genes)


## CORR
t_data = data.frame(t(matrix_out))
corr_mat = cor(t_data)
corrplot::corrplot(corr_mat)

t_data %>% cor_test(vars= t_data$X1556133_s_at , vars2 = c(t_data$X1554112_a_at, t_data$X1569827_at, t_data$X202512_s_at))

test = lm(formula = "TEX36-AS1" ~ genes, data= t_data)

e