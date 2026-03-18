library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(enrichplot)

project_path <- "/home/vaibhav11/rna1/"
in_path <- file.path(project_path)
out_path <- file.path(project_path, "GSEA_res/LLY283")

# Create output directory if it doesn't exist
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Load MSigDB gene sets (C2: Curated gene sets, CGP subcollection)
gene_sets_df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
cat(paste("Loaded", nrow(gene_sets_df), "gene set rows\n"))

# Prepare gene sets as list of gene vectors by pathway name
gene_sets <- gene_sets_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Load your differential expression result file
df <- read.delim(file.path(in_path, 'LLY_vs_DMSO_gene_TE_analysis.txt'), row.names = 1)

# Replace zero or NA p-values to avoid NaNs in ranking
df$pvalue[is.na(df$pvalue)] <- 1
df$pvalue[df$pvalue == 0] <- 1e-300

# Calculate rankings (signed log-pvalue)
rankings <- sign(df$log2FoldChange) * (-log10(df$pvalue))
names(rankings) <- rownames(df)

# Filter rankings to only genes present in the gene sets
genes_in_sets <- unique(unlist(gene_sets))
filtered_genes <- intersect(names(rankings), genes_in_sets)
rankings <- rankings[filtered_genes]

cat(paste("Number of genes in rankings after filtering:", length(rankings), "\n"))
cat("Running fgsea with", length(rankings), "genes\n")

# Sort rankings decreasing
rankings <- sort(rankings, decreasing = TRUE)

# Run fgsea
gsea_res <- fgsea(pathways = gene_sets,
                  stats = rankings,
                  scoreType = 'std',
                  maxSize = 500,
                  minSize = 15)

# Check if fgsea returned results
if (nrow(gsea_res) == 0) {
  cat("fgsea returned no results.\n")
} else {
  print(head(gsea_res))
  cat("Number of significant pathways (padj < 0.05): ", sum(gsea_res$padj < 0.05, na.rm = TRUE), "\n")
}

# If no significant pathways, try a relaxed filter and exploratory plotting
if (sum(gsea_res$padj < 0.05, na.rm = TRUE) == 0) {
  cat("No significant pathways found with padj < 0.05. Trying relaxed threshold and exploratory plots...\n")
  
  # Relaxed filter (padj < 0.25)
  relaxed_pathways <- gsea_res %>% filter(padj < 0.25)
  cat("Number of pathways with padj < 0.25: ", nrow(relaxed_pathways), "\n")
  
  # Pick top 10 pathways by absolute NES regardless of padj
  topPathways <- gsea_res %>% arrange(desc(abs(NES))) %>% head(10) %>% pull(pathway)
  
  # Exploratory enrichment plots
  pdf(file.path(out_path, "top10_pathways_no_sig_filter.pdf"), width = 8, height = 6)
  for (pathway in topPathways) {
    if (!is.null(gene_sets[[pathway]])) {
      p <- plotEnrichment(gene_sets[[pathway]], rankings) + labs(title = pathway)
      print(p)
    }
  }
  dev.off()
}

# Continue with collapsing leadingEdge and saving results as before
if (nrow(gsea_res) > 0) {
  gsea_res$leadingEdge <- sapply(gsea_res$leadingEdge, paste, collapse = ", ")
  write.table(gsea_res, file = file.path(out_path, 'gsea_H.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("No fgsea results to save.\n")
}

# Select top enriched pathways for up and down
