# Load libraries
library(DESeq2)
library(readr)
library(dplyr)
library(tibble)

# Read Telescope counts
counts_MCF10A <- read_csv("MCF10A.csv") %>%
  select(transcript, final_count) %>%
  column_to_rownames(var="transcript")

counts_MDAMB <- read_csv("MDAMB.csv") %>%
  select(transcript, final_count) %>%
  column_to_rownames(var="transcript")

# Rename
colnames(counts_MCF10A) <- "MCF10A"
colnames(counts_MDAMB) <- "MDAMB"

# Union of features
all_features <- union(rownames(counts_MCF10A), rownames(counts_MDAMB))
counts_MCF10A_full <- counts_MCF10A[all_features,, drop=FALSE]
counts_MDAMB_full  <- counts_MDAMB[all_features,, drop=FALSE]
counts_MCF10A_full[is.na(counts_MCF10A_full)] <- 0
counts_MDAMB_full[is.na(counts_MDAMB_full)] <- 0

# Combine matrix
combined_counts <- cbind(counts_MCF10A_full, counts_MDAMB_full)

# Metadata
sample_info <- data.frame(
  row.names = colnames(combined_counts),
  cellType  = factor(c("MCF10A", "MDAMB"))
)

# Create DESeqDataSet (no ~cellType possible with 1 vs 1 → use ~1)
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(combined_counts),
  colData   = sample_info,
  design    = ~ 1
)

# Library size normalization (as in paper)
lib_sizes <- colSums(combined_counts)
geo_mean <- exp(mean(log(lib_sizes)))
sizeFactors(dds) <- lib_sizes / geo_mean

# Normalized counts
norm_counts <- counts(dds, normalized=TRUE)

# Compute CPM cutoff = 0.5
cpm <- sweep(combined_counts, 2, lib_sizes/1e6, "/")
keep <- rowSums(cpm > 0.5) > 0
norm_counts <- norm_counts[keep, ]

# Compute log2FC (MDAMB vs MCF10A)
log2FC <- log2((norm_counts[, "MDAMB"] + 1) / (norm_counts[, "MCF10A"] + 1))

# Build results
res <- data.frame(
  transcript = rownames(norm_counts),
  MCF10A_norm = norm_counts[, "MCF10A"],
  MDAMB_norm  = norm_counts[, "MDAMB"],
  log2FoldChange = log2FC
)

# Filter by |log2FC| > 1 (paper used this with padj < 0.1)
res_filtered <- subset(res, abs(log2FoldChange) > 1)

# Save results
write.csv(res, "MCF10A_vs_MDAMB_all.csv", row.names=FALSE)
write.csv(res_filtered, "MCF10A_vs_MDAMB_filtered.csv", row.names=FALSE)

# Preview
head(res_filtered)

