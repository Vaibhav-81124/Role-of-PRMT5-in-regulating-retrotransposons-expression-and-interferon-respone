library(readr)
library(dplyr)
library(ggplot2)
library(edgeR)
library(patchwork)

# Load merged transcript count data
data <- read_tsv("merged_transcript_counts_all_replicates.txt")

# Prepare count matrix
counts <- as.matrix(data[,-1])
rownames(counts) <- data$transcript

# Assign groups for replicates
group <- c("GSK", "GSK", "LLY", "LLY")

# Calculate CPM
lib_sizes <- colSums(counts)
cpm <- sweep(counts, 2, lib_sizes, FUN = "/") * 1e6

# Calculate mean CPM per treatment
mean_cpm <- data.frame(
  transcript = rownames(counts),
  CPM_GSK = rowMeans(cpm[, group == "GSK", drop = FALSE]),
  CPM_LLY = rowMeans(cpm[, group == "LLY", drop = FALSE])
)

# Top 10 loci by mean CPM for each group
top10_GSK <- mean_cpm %>% arrange(desc(CPM_GSK)) %>% head(10)
top10_LLY <- mean_cpm %>% arrange(desc(CPM_LLY)) %>% head(10)

# GSK barplot
p1 <- ggplot(top10_GSK, aes(x = reorder(transcript, CPM_GSK), y = CPM_GSK)) +
  geom_bar(stat = "identity", fill = "mediumpurple4") +
  coord_flip() +
  labs(title = "Top 10 Highly Expressed Loci - GSK591",
       x = "Locus", y = "Mean CPM") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title = element_text(size = 18))  # <- Increased title size

# LLY barplot
p2 <- ggplot(top10_LLY, aes(x = reorder(transcript, CPM_LLY), y = CPM_LLY)) +
  geom_bar(stat = "identity", fill = "darkseagreen") +
  coord_flip() +
  labs(title = "Top 10 Highly Expressed Loci - LLY283",
       x = "Locus", y = "Mean CPM") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title = element_text(size = 18))  # <- Increased title size


# Combine plots using patchwork
combined_plot <- p1 + p2 + plot_layout(ncol = 2)

# Save combined plot to PDF
ggsave("Top10_Loci_Comparison_GSK_vs_LLY.pdf", combined_plot, width = 25, height = 10, dpi = 600)
ggsave("Top10_Loci_Comparison_GSK_vs_LLY.png", combined_plot, width = 25, height = 10, dpi = 600)

cat("Combined top 10 loci plots saved to 'Top10_Loci_Comparison_GSK_vs_LLY.pdf'.\n")
