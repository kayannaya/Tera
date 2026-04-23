# Install systemPipeR package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("systemPipeR")

# Load systemPipeR package
library(systemPipeR)

# Example command using systemPipeR's system command
system("cutadapt --help")
system("/home/bacteraia/Documents/cutadapt/ --help")
system("/home/bacteraia/Documents/cutadapt/cutadapt --help")

### CUTADAPT
# Load necessary libraries
library(rvest)

# Specify the path to your FastQC HTML report
fastqc_html <- "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/Ctrl_0_1_fastqc.html"

# Read the HTML file
fastqc <- read_html(fastqc_html)

# Extract adapter sequences detected by FastQC using regular expressions
adapter_matches <- fastqc %>%
  html_nodes("div.module_content pre:contains('Adapter')") %>%
  html_text() %>%
  grep("Sequence:", ., value = TRUE)

# Extract adapter sequences using regex
adapter_sequences <- regmatches(adapter_matches, regexpr("Sequence: (.*)", adapter_matches))

# Remove "Sequence: " prefix
adapter_sequences <- gsub("Sequence: ", "", adapter_sequences)

# Print the detected adapter sequences
cat("Detected Adapter Sequences:\n")
print(adapter_sequences)

# Example: List of paths to FastQC HTML files
fastqc_html_files <- c(
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/Ctrl_0_2_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/Ctrl_1_1_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/Ctrl_1_2_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/Sg2_2_1_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/Sg2_2_2_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/SG2_02_1_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/SG2_02_2_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/SG2_8_1_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/SG2_8_2_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/SG2_9_1_fastqc.html",
  "/home/bacteraia/Documents/RNASeq/combined/01.Raw_read/output/SG2_9_1_fastqc.html"
# Add more files as needed
)
fastqc <- read_html(fastqc_html)

adapter_matches <- fastqc %>%
  html_nodes("div.module_content pre:contains('Adapter')") %>%
  html_text() %>%
  grep("Sequence:", ., value = TRUE)
# Extract adapter sequences using regex
adapter_sequences <- regmatches(adapter_matches, regexpr("Sequence: (.*)", adapter_matches))

# Remove "Sequence: " prefix
adapter_sequences <- gsub("Sequence: ", "", adapter_sequences)

# Print the detected adapter sequences
cat("Detected Adapter Sequences:\n")
print(adapter_sequences)

#Reference genome
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")
library(Rsubread)

buildindex(basename="reference_index", reference="/home/bacteraia/Downloads/GCA_000001405.29_GRCh38.p14_genomic.fasta")

#Genome directory is built in linux termnial using the code below
#STAR --runMode genomeGenerate\
#--genomeDir 05.STAR/ref/\
#--genomeFastaFiles 05.STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa\
#--sjdbGTFfile 05.STAR/Homo_sapiens.GRCh38.112.gtf\
#--runThreadN 16

#alignment, done in parallel to conserve time. However, will be system intensive.
library(BiocParallel)
library(STAR)
gen_ref <- "./05.STAR/ref"
STAR_res <- "./05.STAR/res"

samples <- c("Ctrl_0", "Ctrl_1", "Sg2_2", "Sg2_02", "Sg2_8", "Sg2_9")
fastq_files <- list(
  c("./02.Clean_read/Ctrl_0_trimmed_1.fq", ",/02.Clean_read/Ctrl_0_trimmed_2.fq"),
  c("./02.Clean_read/Ctrl_1_trimmed_1.fq", ",/02.Clean_read/Ctrl_1_trimmed_2.fq"),
  c("./02.Clean_read/Sg2_2_trimmed_1.fq", ",/02.Clean_read/Sg2_2_trimmed_2.fq"),
  c("./02.Clean_read/Sg2_02_trimmed_1.fq", ",/02.Clean_read/Sg2_02_trimmed_2.fq"),
  c("./02.Clean_read/Sg2_8_trimmed_1.fq", ",/02.Clean_read/Sg2_8_trimmed_2.fq"),
  c("./02.Clean_read/Sg2_9_trimmed_1.fq", ",/02.Clean_read/Sg2_9_trimmed_2.fq")
)

# Function to run STAR for each sample
run_star <- function(sample_name, fastq_pair) {
  star_cmd <- paste("STAR",
                    "--runThreadN 3",  # Use 3 threads per sample
                    "--genomeDir", gen_ref,
                    "--readFilesIn", fastq_pair[1], fastq_pair[2], 
                    "--readFilesCommand zcat", 
                    "--outFileNamePrefix", file.path(STAR_res, sample_name),
                    "--outSAMtype BAM SortedByCoordinate")
  system(star_cmd)
}

# Run STAR in parallel using mclapply()
mclapply(1:length(samples), function(i) {
  run_star(samples[i], fastq_files[[i]])
}, mc.cores = 10)  # Using 6 cores in parallel
