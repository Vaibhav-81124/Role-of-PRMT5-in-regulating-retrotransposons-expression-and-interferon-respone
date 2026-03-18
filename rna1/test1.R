library(biomaRt)
library(fgsea)
library(msigdbr)

# Step 1: Set up biomaRt connection (using Ensembl's gene annotation)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Your gene list (both RefSeq and Ensembl IDs)
gene_ids <- c(NM_001034077.4,NM_001042422.2,NM_001042423.2,NM_001085.4,NM_001199579.1,NM_001199580.1,NM_001206950.1,NM_001206951.1,NM_001206952.1,NM_001255.2,NM_001258311.1,NM_002229.2,NM_003519.3,NM_003524.2,NM_003538.3,NM_003539.3,NM_003544.2,NM_004207.3,NM_005030.5,NM_005764.3,NM_006427.3,NM_014501.2,NM_015374.2,NM_018942.2,NM_020650.2,NM_021066.2,NM_021709.2,NM_031454.1,NM_133467.2,NM_145030.2,NR_002995.1,NR_027033.2,NR_110479.1,NR_119376.1,NR_132278.1,XM_005259089.3,XM_005268325.3,XM_006716616.2,XM_006718223.2,XM_006720869.2,XM_011515535.2,XM_011523607.1,XM_011526752.1,XM_011527143.1,XM_011530104.2,XM_011530105.2,XM_011531191.2,XM_011536360.2,XM_017000537.1,XM_017002330.1,XM_017012655.1,XM_017027023.1,XR_929661.2)

# Split the IDs into two categories: RefSeq (NM, NR) and Ensembl (XM, XR)
refseq_ids <- gene_ids[grepl("^NM_|^NR_", gene_ids)]  # RefSeq IDs
ensembl_ids <- gene_ids[grepl("^XM_|^XR_", gene_ids)]  # Ensembl IDs

# Step 2: Map RefSeq IDs to Entrez Gene IDs using biomaRt
refseq_mapping <- getBM(
  attributes = c("refseq_mrna", "entrezgene_id", "hgnc_symbol"),
  filters = "refseq_mrna", 
  values = refseq_ids,
  mart = ensembl
)

# Map Ensembl IDs to Entrez Gene IDs using biomaRt
ensembl_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "entrezgene_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id", 
  values = ensembl_ids,
  mart = ensembl
)

# Combine the results
combined_mapping <- rbind(refseq_mapping, ensembl_mapping)

# View the mapping result (transcript IDs to Entrez Gene IDs)
print(combined_mapping)

# Step 3: Prepare Gene Expression Data (if you have it)
# Example data with log fold changes (assuming this is your expression data)
gene_expression <- data.frame(
  gene = c("NM_001034077.4", "XM_005259089.3", "NR_002995.1"),
  logFC = c(2.1, -1.5, 3.0)
)

# Map the gene expression data to Entrez Gene IDs
gene_expression$entrezgene_id <- combined_mapping$entrezgene_id[match(gene_expression$gene, combined_mapping$refseq_mrna)]

# Ensure no missing Entrez Gene IDs
gene_expression <- gene_expression[!is.na(gene_expression$entrezgene_id), ]

# Sort by logFC to create the ranked list
ranked_genes <- gene_expression[order(gene_expression$logFC, decreasing = TRUE), ]

# Create a named vector for fgsea
ranked_genes <- setNames(ranked_genes$logFC, ranked_genes$entrezgene_id)

# Step 4: Load Gene Sets from MSigDB
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- hallmark_gene_sets %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarize(genes = list(unique(gene_symbol))) %>%
  dplyr::pull(genes, gs_name)

# Step 5: Run fgsea
fgsea_results <- fgsea(pathways = gene_sets, 
                       stats = ranked_genes, 
                       minSize = 15, 
                       maxSize = 500, 
                       nPermSimple = 1000)

# View the results
head(fgsea_results)

# Optional: Plot the results (Top 10 Enriched Pathways)
library(ggplot2)
ggplot(fgsea_results[order(pval), ][1:10, ], aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Enriched Pathways", y = "Normalized Enrichment Score (NES)", x = "Gene Set")
