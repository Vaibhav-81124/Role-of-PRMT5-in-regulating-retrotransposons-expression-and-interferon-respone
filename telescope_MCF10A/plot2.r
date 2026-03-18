library(tidyverse)
library(patchwork)

# Load CPM normalized counts file
df <- read.csv("cpm_normalized_counts.csv")

# Select top 10 loci by CPM for each cell type
top_mcfa <- df %>% arrange(desc(MCF10A)) %>% slice(1:10)
top_mdamb <- df %>% arrange(desc(MDAMB436)) %>% slice(1:10)

# Plot function (no error bars)
plot_top10_single <- function(stats_df, cell_type, col_name) {
  ggplot(stats_df, aes(x = reorder(transcript, !!sym(col_name)), y = !!sym(col_name))) +
    geom_bar(stat = "identity", fill = "darkslateblue") +
    coord_flip() +
    labs(title = paste("Top 10 Highly Expressed Loci -", cell_type),
         x = "Transcript",
         y = "CPM") +
    theme_minimal()
}

# Draw plots
p1 <- plot_top10_single(top_mcfa, "MCF10A", "MCF10A")
p2 <- plot_top10_single(top_mdamb, "MDAMB436", "MDAMB436")

# Combine with patchwork
combined_plot <- p1 | p2   # side-by-side layout
# Or use p1 / p2 for stacked vertically

# Save combined plot
ggsave("Top10_Loci_Combined.png", combined_plot, width = 14, height = 6)

