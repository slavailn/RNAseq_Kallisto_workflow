library(tximport)
library(DESeq2)
library(ensembldb)
library(rtracklayer)

setwd(<PROJECT_DIR>)
list.files("kallisto_quant/")

# Create transcript to gene map using the original gene annotation file
# in GTF format. Here we are using the same gtf as one in used to create
# transctipt transcript sequences
gtf <- rtracklayer::import("docs/genes.gtf")
tx2gene <- unique(data.frame(
  TXNAME = mcols(gtf)$transcript_id,
  GENEID = mcols(gtf)$gene_id,
  stringsAsFactors = FALSE
))
head(tx2gene)

# Import kallisto abundance estimates
# Change pattern recognition in list files and sample names as appropriate for 
# your projects
files <- list.files("kallisto_quant", pattern = "abundance.tsv$", 
                    recursive = TRUE, full.names = TRUE)
dirs <- basename(dirname(files))
print(dirs)
samples <- sub(".*\\.", "", dirs)
# Note, that files must be a named vector for tximport object to have
# column names. Here, I gave files vector names corresponding to sample
# names in my project.
names(files) <- samples

txi <- tximport(files, type = "kallisto", txIn = TRUE,
                txOut = FALSE, countsFromAbundance = "lengthScaledTPM",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)
head(txi$abundance)

# Create DESeq2 object from tximport output
# First, load metadata
metadata <- read.csv("docs/metadata.csv", header = T)
head(metadata)
rownames(metadata) <- metadata$SampleID

# Re-order rows in metadata to match columns sample order in txi object
meta_ordered <- 
  metadata[match(colnames(txi$abundance), rownames(metadata)), , drop = FALSE]
head(meta_ordered)
identical(rownames(meta_ordered), colnames(txi$abundance))
print(data.frame(rownames(meta_ordered), 
                 colnames(txi$abundance)))

# Create DESeq2 object, save intermediate files
dds <- DESeqDataSetFromTximport(txi, colData = meta_ordered, 
                                design = ~Treatment)
save(dds, file = "DESeqObject.RData")

# Run normalization and variance stabilization
vsd <- varianceStabilizingTransformation(dds)
head(assay(vsd))

# ------ END ------ #
# Kallisto output is now loaded into R and turned into DESeq2 object.
# DESeq2 object was normalized and variance stabilized. This type of data
# is suitable for quality control, clusering and machine learning.





save(vsd, file = "DESeq2_VST.RData")
