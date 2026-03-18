# Load libraries
library(ggplot2)
library(dplyr)
library(patchwork)   # for combining plots

# Read DEG table
deg_df <- read.csv("MCF10A_vs_MDAMB_filtered1.csv")

# Histogram of log2 fold changes
p1 <- ggplot(deg_df, aes(x = log2FoldChange)) +
  geom_histogram(bins = 60, fill = "darkviolet", alpha = 0.7) +
  theme_minimal() +
  xlab("Log2 Fold Change (MCF10A vs MDAMB)") +
  ggtitle("Histogram of Log2 Fold Changes")

# Add average expression for MA-like plot
deg_df <- deg_df %>%
  mutate(avgExp = (MCF10A_norm + MDAMB_norm) / 2)

p2 <- ggplot(deg_df, aes(x = log10(avgExp + 1), y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  xlab("Log10 Mean Normalized Expression") +
  ylab("Log2 Fold Change") +
  ggtitle("MA Plot-like for LINE1 Expression")

# Top 10 upregulated transcripts
top_up <- deg_df %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10) %>%
  mutate(transcript = reorder(transcript, log2FoldChange))

p3 <- ggplot(top_up, aes(x = transcript, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "coral2") +
  coord_flip() +
  theme_minimal() +
  xlab("LINE1 Transcript") +
  ylab("Log2 Fold Change") +
  ggtitle("Top 10 Upregulated LINE1 locus (MDAMB vs MCF10A)")

# Top 10 downregulated transcripts
top_down <- deg_df %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10) %>%
  mutate(transcript = reorder(transcript, log2FoldChange))

p4 <- ggplot(top_down, aes(x = transcript, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "cornflowerblue") +
  coord_flip() +
  theme_minimal() +
  xlab("LINE1 Transcript") +
  ylab("Log2 Fold Change") +
  ggtitle("Top 10 Downregulated LINE1 locus (MCF10A vs MDAMB)")

# Combine all plots into 2x2 grid
combined_plot <- (p1 | p2) / (p3 | p4)

# Save all plots into a single PDF
ggsave("LINE1_DE_Plots.pdf", combined_plot, width = 14, height = 10)

# If you also want a PNG:
ggsave("LINE1_DE_Plots.png", combined_plot, width = 14, height = 10, dpi = 300)
