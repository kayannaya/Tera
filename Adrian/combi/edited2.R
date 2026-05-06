#Deseq2
BiocManager::install("DESeq2", force = TRUE)
library(DESeq2)

counts <- read.table("./06_DESeq2/counts_2.txt", header = TRUE, row.names = 1, sep = "\t")
head(counts)

sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c("control", "control", "treatment", "treatment", "treatment")
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ condition)

dds <- DESeq(dds)
head(dds)
res <- results(dds)
head(res)
write.csv(as.data.frame(res), file = "./06_DESeq2/deseq2_results_unfiltered_2.csv")

res_filtered <- res[which(res$padj < 0.05), ]
head(res_filtered)
write.csv(as.data.frame(res_filtered), file = "./06_DESeq2/deseq2_results_filtered_2.csv")

#Ensmbl to Gene ID
BiocManager::install(update = TRUE, ask = FALSE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("AnnotationDbi", force = TRUE)
library(org.Hs.eg.db)
res.df <- as.data.frame(res)
res.df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res.df),
                        keytype = "ENSEMBL",
                        column = "SYMBOL",
                        multiVals = "first")
head(res.df)
write.csv(as.data.frame(res.df), file = "./06_DESeq2/deseq2_results_unfiltered_geneName.csv")

res_f.df <- as.data.frame(res_filtered)
res_f.df$symbol <- mapIds(org.Hs.eg.db,
                          keys = rownames(res_f.df),
                          keytype = "ENSEMBL",
                          column = "SYMBOL",
                          multiVals = "first")
head(res_f.df)
write.csv(as.data.frame(res_f.df), file = "./06_DESeq2/deseq2_results_filtered_geneName.csv")

#Plot
BiocManager::install("ggplot2", force = TRUE)
library(ggplot2)

gene <- "ENSG00000171889"
normalized_counts <- counts(dds, normalized = TRUE) [gene, colnames(dds) ]
plot_data <- data.frame(
  condition = colData(dds)$condition,  
  counts = normalized_counts
)

ggplot(plot_data, aes(x = condition, y = counts)) +
  geom_boxplot() +
  labs(title = paste("Expression of", gene),
       x = "Condition",
       y = "Normalized Counts") +
  theme_minimal()
 #Violin
ggplot(plot_data, aes(x = condition, y = counts)) +
  geom_violin() +
  labs(title = paste("Expression of MIR31HG", gene),
       x = "Condition",
       y = "Normalized Counts") +
  theme_minimal()
#dot plot
ggplot(plot_data, aes(x = condition, y = counts)) +
  geom_jitter(width = 0.2) +
  labs(title = paste("Expression of", gene),
       x = "Condition",
       y = "Normalized Counts") +
  theme_minimal()
