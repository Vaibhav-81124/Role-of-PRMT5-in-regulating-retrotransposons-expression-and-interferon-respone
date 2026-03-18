# Load necessary libraries
library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(enrichplot)
library(patchwork)  # Added patchwork for combining plots

# Define project paths
project_path <- "/home/vaibhav11/GSEA1/GSK591_vs_DMSO_1d"
in_path <- file.path(project_path)
out_path <- file.path(project_path, "GSEA_results/H")  

# msigdb data input
gene_sets_df <- msigdbr(species = "Homo sapiens", collection = "H")
cat(paste("loaded", length(gene_sets_df), "gene sets\n"))

# Split gene sets
gene_sets <- gene_sets_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Read differential expression file
df <- read.delim(file.path(in_path, 'diff-exp-fdr.txt'), row.names = 1)

# Create rankings for GSEA
rankings <- sign(df$log2FoldChange) * (-log10(df$pvalue))
names(rankings) <- rownames(df)
rankings <- sort(rankings, decreasing = TRUE)

# Perform GSEA
gsea_res <- fgsea(pathways = gene_sets,
                  stats = rankings,
                  scoreType = 'std',
                  maxSize = 500,
                  minSize = 15)
gsea_res$leadingEdge <- sapply(gsea_res$leadingEdge, paste, collapse = ", ")
write.table(gsea_res, file = file.path(out_path, 'gsea_H.txt'), sep = "\t", quote = F, row.names = F)

# Upregulated and Downregulated pathways
number_of_top_pathways_up <- 10
number_of_top_pathways_DN <- 10
toppathwaysUP <- gsea_res[NES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]
toppathwaysDN <- gsea_res[NES < 0][head(order(padj), n = number_of_top_pathways_DN), pathway]

topPathways <- c(toppathwaysUP, rev(toppathwaysDN))

# Bar Plot for Top Upregulated Pathways
bar_data <- gsea_res[gsea_res$pathway %in% toppathwaysUP, ] %>%
  arrange(desc(NES)) %>%
  mutate(pathway = factor(pathway, levels = rev(pathway)))

bar_plot <- ggplot(bar_data, aes(x = NES, y = pathway, fill = NES)) +
  geom_col() +
  scale_fill_gradient2(low = "white", high = "blue", midpoint = 0) +
  theme_minimal(base_size = 14) +
  labs(title = "Top UP-Regulated Enriched Pathways (Bar Plot)",
       x = "Normalized Enrichment Score (NES)",
       y = "Pathway")

# Dot Plot for Top Upregulated Pathways
dot_data <- gsea_res[gsea_res$pathway %in% toppathwaysUP, ] %>%
  arrange(padj) %>%
  mutate(pathway = factor(pathway, levels = rev(pathway)))

dot_plot <- ggplot(dot_data, aes(x = NES, y = pathway, size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "white", high = "blue", midpoint = 0) +
  theme_minimal(base_size = 14) +
  labs(title = "Top UP-Regulated Enriched Pathways (Dot Plot)",
       x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       size = "-log10(adj. p-value)")

# Enrichment Plots for Upregulated Pathways (first one as an example)
up_enrichment_plots <- list()
for (path in toppathwaysUP) {
  if (!is.null(gene_sets[[path]])) {
    p <- plotEnrichment(gene_sets[[path]], rankings) + labs(title = path)
    up_enrichment_plots[[path]] <- p
  }
}

# Select first enrichment plot for UP-regulated pathways
first_up_enrichment_plot <- up_enrichment_plots[[1]]

# Enrichment Plots for Downregulated Pathways (first one as an example)
down_enrichment_plots <- list()
for (path in toppathwaysDN) {
  if (!is.null(gene_sets[[path]])) {
    p <- plotEnrichment(gene_sets[[path]], rankings) + labs(title = path)
    down_enrichment_plots[[path]] <- p
  }
}

# Select first enrichment plot for Down-regulated pathways
first_down_enrichment_plot <- down_enrichment_plots[[1]]

# Pathway plot (GSEA Table)
pathway_plot <- plotGseaTable(gene_sets[topPathways], stats = rankings, fgseaRes = gsea_res, gseaParam = 1)

# Combine all plots using patchwork
combined_plot <- (bar_plot | dot_plot) / (first_up_enrichment_plot | first_down_enrichment_plot) / pathway_plot

# Save the combined plot as a PDF
ggsave(file = file.path(out_path, "combined_plot.pdf"), plot = combined_plot, width = 20, height = 15)

# You can also display it in RStudio viewer
combined_plot
