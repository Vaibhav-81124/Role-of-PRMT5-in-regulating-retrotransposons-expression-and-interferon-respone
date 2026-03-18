library(readr)
library(DESeq2)
library(edgeR)

# Load merged transcript count data
merged_counts <- read_tsv("merged_transcript_counts_all_replicates.txt")

# Extract count matrix and set transcript names as rownames
count_matrix <- as.matrix(merged_counts[,-1])
rownames(count_matrix) <- merged_counts$transcript

# Calculate library sizes
lib_sizes <- colSums(count_matrix)

# CPM normalization using edgeR
cpm_counts <- cpm(count_matrix, lib.size = lib_sizes)

# Create a sample metadata dataframe assuming two treatments with two replicates each
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  Treatment = factor(c("GSK_rep1", "GSK_rep2", "LLY_rep1", "LLY_rep2"))
)

# Create DESeq2 dataset for VST normalization
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Treatment)

# Prefilter loci with low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2 to estimate size factors before VST
dds <- estimateSizeFactors(dds)

# VST normalization
vst_data <- varianceStabilizingTransformation(dds, blind = TRUE)

# Extract VST assay matrix
vst_mat <- assay(vst_data)

# Save normalized outputs to files
write.csv(cpm_counts, file = "CPM_normalized_counts.csv")
write.csv(vst_mat, file = "VST_normalized_counts.csv")

cat("CPM and VST normalization complete. Outputs saved as CSV files.\n")
