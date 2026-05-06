#GO
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), force = TRUE)
BiocManager::install(c("readr", "dplyr", "AnnotationDbi"), force = TRUE)
##Installing DOSE
BiocManager::install("DOSE", force = TRUE)
remove.packages("DOSE")

## Installing org hs eg db
BiocManager::install("org.Hs.eg.db", force = TRUE)

## Installing the clusterProfiler
BiocManager::install("clusterProfiler", force = TRUE)

##Installing the Pathview
BiocManager::install("pathview", force = TRUE)

library(clusterProfiler)
remove.packages("clusterProfiler")
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(pathview)
keytypes(org.Hs.eg.db)

deseq2_results <- read.csv("./06_DESeq2/deseq2_results_unfiltered_2.csv")
write.table(write.table(deseq2_results, "./06_DESeq2/GO.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE))
genes <- read.table("./06_DESeq2/GO.txt", header = TRUE, sep = "\t")

str(genes)

# Extract Ensembl IDs
ensembl_ids <- genes$ENSEMBL

# Ensure ensembl_ids is a character vector
ensembl_ids <- as.character(ensembl_ids)

# Convert Ensembl IDs to Entrez IDs
entrez_ids <- bitr(ensembl_ids, 
                   fromType = 'ENSEMBL', 
                   toType = 'ENTREZID', 
                   OrgDb = org.Hs.eg.db)

# View the conversion result
head(entrez_ids)

gene.df <- entrez_ids$ENTREZID
head(gene.df)
####GO classification
library(clusterProfiler)
head(genes)
go <- enrichGO(gene = genes, 
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "ENSEMBL",
                       ont = "BP",  # Biological Process; use "CC" for Cellular Component, "MF" for Molecular Function
                       qvalueCutoff = 0.05)
