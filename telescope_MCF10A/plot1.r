library(tidyverse)
library(patchwork)

# Load normalized CPM counts
df <- read.csv("cpm_normalized_counts.csv")

# If family column missing, add manually based on transcript name
df <- df %>%
  mutate(family = ifelse(grepl("L1FLnI", transcript), "L1FLnI", "L1ORF2"))

# Convert from wide to long format
df_long <- df %>%
  pivot_longer(
    cols = -c(transcript, family),
    names_to = "sample",
    values_to = "CPM"
  )

# --- Family-wise proportion per sample ---
p1 <- df_long %>%
  group_by(sample, family) %>%
  summarise(totalCPM = sum(CPM), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(propCPM = totalCPM / sum(totalCPM)) %>%
  ggplot(aes(x = sample, y = propCPM, fill = family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  ylab("Proportion of LINE-1 Expression") +
  xlab("Sample") +
  ggtitle("Proportion of LINE-1 Families per Sample")

# --- Boxplot: Number of expressed LINE1 elements per family ---
p2 <- df_long %>%
  filter(CPM > 0.5) %>%
  group_by(sample, family) %>%
  summarise(n_expressed = n(), .groups = "drop") %>%
  ggplot(aes(x = family, y = n_expressed, fill = family)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  theme_minimal() +
  ylab("Number of Expressed LINE-1 Elements") +
  xlab("LINE-1 Family") +
  ggtitle("Expressed LINE-1 Elements per Family Across Samples")

# --- Combine with patchwork ---
combined_plot <- p1 / p2  # stack vertically (use p1 | p2 for side-by-side)

# Save to PNG
ggsave("LINE1_family_expression_plots.png", combined_plot,
       width = 10, height = 10, dpi = 600)
