library(readr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(pvclust)

# Load raw count files (no header assumed)
mcf10a <- read_tsv("MCF10A-telescope_report.tsv", col_names = FALSE, show_col_types = FALSE)
mdamb436 <- read_tsv("MDAMB436-telescope_report.tsv", col_names = FALSE, show_col_types = FALSE)

# Select locus and raw counts
mcf10a_counts <- mcf10a %>% select(locus = X1, count_MCF10A = X3)
mdamb436_counts <- mdamb436 %>% select(locus = X1, count_MDAMB436 = X3)

# Merge counts on locus
combined_counts <- full_join(mcf10a_counts, mdamb436_counts, by = "locus")

# Convert count columns to numeric then replace NA with 0
combined_counts <- combined_counts %>%
  mutate(across(starts_with("count"), as.numeric)) %>%
  mutate(across(starts_with("count"), ~replace_na(., 0)))

# Save filtered raw counts
write.csv(combined_counts, "filtered_raw_counts.csv", row.names = FALSE)

# Calculate library sizes and CPM
library_sizes <- c(sum(combined_counts$count_MCF10A), sum(combined_counts$count_MDAMB436))
cpm_MCF10A <- (combined_counts$count_MCF10A / library_sizes[1]) * 1e6
cpm_MDAMB436 <- (combined_counts$count_MDAMB436 / library_sizes[2]) * 1e6
expressed_loci <- (cpm_MCF10A >= 0.5) | (cpm_MDAMB436 >= 0.5)
filtered_counts <- combined_counts[expressed_loci, ]

# Save filtered counts
write.csv(filtered_counts, "expressed_filtered_counts.csv", row.names = FALSE)

# Create count matrix for DESeq2
count_matrix <- filtered_counts %>%
  select(starts_with("count")) %>%
  as.matrix()
rownames(count_matrix) <- filtered_counts$locus

# Sample metadata with intercept-only design (no replicates)
col_data <- data.frame(row.names = colnames(count_matrix))

# DESeq2 dataset with design ~1 (intercept only)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ 1)

# Variance stabilizing transformation in blind mode
vst_data <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_mat <- assay(vst_data)

# Remove rows with any NAs to avoid pvclust errors
vst_mat_clean <- vst_mat[complete.cases(vst_mat), ]

# Save VST normalized counts
write.csv(as.data.frame(vst_mat_clean), "vst_normalized_counts_clean.csv")

# Hierarchical clustering with bootstrap resampling
pc <- pvclust(vst_mat_clean, method.hclust = "average", method.dist = "correlation", nboot = 1000)

# Save clustering result object
save(pc, file = "pvclust_results.RData")

# Plot dendrogram with bootstrap confidence values
plot(pc)
pvrect(pc, alpha = 0.95)

cat("VST normalization and clustering complete. Cleaned matrix saved, clustering plotted.\n")

