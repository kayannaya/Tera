# Read the file without the '.txt' extension
data <- read.table("./05_Read_alignment/counts_out.txt", header = TRUE, sep = "\t")

# Check the first few rows
head(data)

data <- data[ , -2:-6]

print(colnames(data))
colnames(data) <- substr(colnames(data), 15, nchar(colnames(data)))
colnames(data) <- substr(colnames(data), 1, 7)

# Save the data frame as a CSV file
write.csv(data, "./06_DESeq2/data.csv", row.names = FALSE)
write.table(data, "./06_DESeq2/counts.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Deseq2
BiocManager::install("DESeq2", force = TRUE)
library(DESeq2)

counts <- read.table("./06_DESeq2/counts.txt", header = TRUE, row.names = 1, sep = "\t")
head(counts)

sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c("control", "control", "treatment", "treatment", "treatment", "treatment")
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ condition)

dds <- DESeq(dds)
head(dds)
res <- results(dds)
head(res)
write.csv(as.data.frame(res), file = "./06_DESeq2/deseq2_results_unfiltered.csv")

res_filtered <- res[which(res$padj < 0.05), ]
head(res_filtered)
write.csv(as.data.frame(res), file = "./06_DESeq2/deseq2_results_filtered.csv")

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
normalized_counts <- counts(dds, normalized = TRUE)
plot_data <- data.frame(
  condition = colData(dds)$condition,  # Assuming you have a 'condition' column in colData
  counts = gene_counts
)

ggplot(plot_data, aes(x = condition, y = counts)) +
  geom_boxplot() +
  labs(title = paste("Expression of MIR31HG", gene),
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
