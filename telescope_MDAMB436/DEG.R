library(readr)
library(DESeq2)
library(dplyr)

# Load merged counts CSV or TSV
merged_counts <- read_tsv("merged_transcript_counts_all_replicates.txt")

# Prepare count matrix (remove transcript column, set rownames)
count_matrix <- as.matrix(merged_counts[,-1])
rownames(count_matrix) <- merged_counts$transcript

# Create sample metadata matching column names exactly
sample_metadata <- data.frame(
  condition = factor(c("GSK", "GSK", "LLY", "LLY")),
  row.names = colnames(count_matrix)
)

# Check names match
stopifnot(all(colnames(count_matrix) == rownames(sample_metadata)))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_metadata,
                              design = ~ condition)

# Prefilter low count loci
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Add normalized counts means per condition for visualization
norm_counts <- counts(dds, normalized=TRUE)
resDF <- as.data.frame(res) %>%
  mutate(transcript = rownames(res),
         LLY_norm = rowMeans(norm_counts[, sample_metadata$condition == "LLY"]),
         GSK_norm = rowMeans(norm_counts[, sample_metadata$condition == "GSK"]))

# Save results for plotting and further analysis
write.csv(resDF, "LLY_vs_GSK_DESeq2_results.csv", row.names = FALSE)

cat("DESeq2 analysis complete. Results saved as 'LLY_vs_GSK_DESeq2_results.csv'.\n")
