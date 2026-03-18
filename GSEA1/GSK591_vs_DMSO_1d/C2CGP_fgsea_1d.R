library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(enrichplot)

project_path <- "/home/vaibhav11/GSEA1/GSK591_vs_DMSO_1d"
in_path <- file.path(project_path)
out_path <- file.path(project_path, "GSEA_results/C2CGP")  
#msigdb data input
gene_sets_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
cat(paste("loaded", length(gene_sets_df), "gene sets\n"))

#df
gene_sets <- gene_sets_df %>% split(x = .$gene_symbol, f = .$gs_name)
df <- read.delim(file.path(in_path, 'diff-exp-fdr.txt'), row.names = 1)

#rank file
rankings <- sign(df$log2FoldChange) * (-log10(df$pvalue))
names(rankings) <- rownames(df)
rankings <- sort(rankings, decreasing = TRUE)
write.table(rankings,file=file.path(out_path,"ranks.txt"),sep = "\t")

#gsea 
gsea_res <- fgsea(pathways = gene_sets,
                  stats = rankings,
                  scoreType = 'std',
                  maxSize = 500,
                  minSize = 15)
gsea_res$leadingEdge <- sapply(gsea_res$leadingEdge, paste, collapse = ", ")
write.table(gsea_res, file = file.path(out_path, 'gsea_C2CGP.txt'), sep = "\t", quote = F, row.names = F)

# upregulated and dowregulated pathways
number_of_top_pathways_up <- 10
number_of_top_pathways_DN <- 10
toppathwaysUP <- gsea_res[NES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]
toppathwaysDN <- gsea_res[NES < 0][head(order(padj), n = number_of_top_pathways_DN), pathway]

topPathways <- c(toppathwaysUP, rev(toppathwaysDN))

# gseaplot
pdf(file = file.path(out_path, paste0("C2CGP_", 'pathway.pdf')), width = 20, height = 15)
plotGseaTable(gene_sets[topPathways], stats = rankings, fgseaRes = gsea_res, gseaParam = 1)
dev.off()

# Enrichment plot for upregulated
toppathwaysUP <- unique(toppathwaysUP)
pdf(file = file.path(out_path, "enrichment_plots_up_pathways.pdf"), width = 8, height = 6)
for (path in toppathwaysUP) {
  if (!is.null(gene_sets[[path]])) {
    p <- plotEnrichment(gene_sets[[path]], rankings) + labs(title = path)
    print(p)
  }
}
dev.off()

# Enrichment plot for downregulated
toppathwaysDN <- unique(toppathwaysDN)
pdf(file = file.path(out_path, "enrichment_plots_DN_pathways.pdf"), width = 8, height = 6)
for (path in toppathwaysDN) {
  if (!is.null(gene_sets[[path]])) {
    p <- plotEnrichment(gene_sets[[path]], rankings) + labs(title = path)
    print(p)
  }
}
dev.off()

# Barplot
bar_data <- gsea_res[gsea_res$pathway %in% toppathwaysUP, ] %>%
  arrange(desc(NES)) %>%
  mutate(pathway = factor(pathway, levels = rev(pathway)))

ggplot(bar_data, aes(x = NES, y = pathway, fill = NES)) +
  geom_col() +
  scale_fill_gradient2(low = "white", high = "blue", midpoint = 0) +
  theme_minimal(base_size = 14) +
  labs(title = "Top UP-Regulated Enriched Pathways (Bar Plot)",
       x = "Normalized Enrichment Score (NES)",
       y = "Pathway")

ggsave(file = file.path(out_path, "barplot_toppathwaysUP.pdf"), width = 15, height = 6)

# Dotplot
dot_data <- gsea_res[gsea_res$pathway %in% toppathwaysUP, ] %>%
  arrange(padj) %>%
  mutate(pathway = factor(pathway, levels = rev(pathway)))

ggplot(dot_data, aes(x = NES, y = pathway, size = -log10(padj), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "white", high = "blue", midpoint = 0) +
  theme_minimal(base_size = 14) +
  labs(title = "Top UP-Regulated Enriched Pathways (Dot Plot)",
       x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       size = "-log10(adj. p-value)")

ggsave(file = file.path(out_path, "dotplot_toppathwaysUP.pdf"), width = 15, height = 6)

# Cnet
up_gene_sets <- gene_sets[toppathwaysUP]
pdf(file = file.path(out_path, "cnetplot_upregulated.pdf"), width = 120, height = 120)
cnetplot(up_gene_sets, foldChange = rankings, showCategory = 5, circular = FALSE, max.overlap = 50000, layout = "kk")
dev.off()

# Cnet
down_gene_sets <- gene_sets[toppathwaysDN]
pdf(file = file.path(out_path, "cnetplot_downregulated.pdf"), width = 60, height = 60)
cnetplot(down_gene_sets, foldChange = rankings, showCategory = 5, circular = FALSE, max.overlap = 50000, layout = "kk")
dev.off()


