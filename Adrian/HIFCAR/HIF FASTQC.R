if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("Rsamtools", force=TRUE) 
BiocManager::install("GenomicAlignments", force=TRUE)
BiocManager::install("GenomicFeatures", force=TRUE)
BiocManager::install("Rsamtools", force=TRUE)
BiocManager::install("DESeq2", force=TRUE)
BiocManager::install("ggplot2", force=TRUE)
BiocManager::install("pheatmap", force=TRUE)
BiocManager::install("cutadapt", force=TRUE)
library (system)
library(FastQC)

###FASTQC
fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/Ctrl_0_1.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/Ctrl_0_2.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/Ctrl_1_1.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/Ctrl_1_2.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/Sg2_2_1.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/Sg2_2_2.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/SG2_02_1.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/SG2_02_2.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/SG2_8_1.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/SG2_8_2.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/SG2_9_1.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

fastq_file <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/SG2_9_2.fq.gz"
output_dir <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output"
fastqc_command <- paste("fastqc", fastq_file, "-o", output_dir)
system(fastqc_command)

