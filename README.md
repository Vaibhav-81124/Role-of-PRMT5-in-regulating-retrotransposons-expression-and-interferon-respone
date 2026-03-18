# PRMT5 and Retrotransposon Regulation

**Investigating the Role of PRMT5 in Regulating Retrotransposon Expression and Interferon Response**  
SRFP 2025 Project – Indian Academy of Sciences  
Author: Vaibhav D  

---

## Overview
This project investigates the role of Protein Arginine Methyltransferase 5 (PRMT5) in regulating transposable elements and interferon response pathways using integrative bioinformatics approaches.  
The study focuses on how PRMT5 inhibition influences:
- Retrotransposon expression
- Epigenetic regulation
- Immune signaling pathways

---

## Background
PRMT5 is a key epigenetic regulator involved in histone methylation, RNA splicing, and immune modulation. Transposable elements (TEs), which constitute nearly half of the human genome, are typically silenced through epigenetic mechanisms.  

Disruption of PRMT5 function can:
- Reactivate transposable elements such as endogenous retroviruses (ERVs) and LINE-1 elements
- Trigger innate immune pathways through cytosolic nucleic acid sensing

---

## Objectives
- Analyze the impact of PRMT5 inhibition on retrotransposon expression
- Evaluate changes in interferon-related pathways
- Identify key transposable element families and loci affected
- Study PRMT5 binding patterns on LINE-1 elements
- Quantify TE expression at both family and locus levels

---

## Methodology

### Gene Set Enrichment Analysis (GSEA)
- Dataset: GEO (GSE236930)
- Tool: fgsea (R)
- Gene sets: Hallmark, C2 (KEGG/Reactome), C5 (GO)
- Focus: interferon and cell-cycle pathways

### RNA-seq Analysis
- Quality control: FastQC
- Trimming: Trim Galore
- Alignment: STAR
- Differential expression: DESeq2

### Transposable Element Analysis
- Tool: TEtranscripts
- Uses expectation-maximization to handle multi-mapped reads
- Outputs: differential TE expression and class distribution

### Locus-Level TE Quantification
- Tool: Telescope
- Resolves multi-mapped reads to specific genomic loci

### ChIP-seq Analysis
- Tools: Bowtie2, Samtools, HOMER
- Objective: identify PRMT5 binding on LINE-1 elements
- Visualization: IGV

### Primer Design
- Designed PRMT5 construct for FLAG-tag expression
- In-Fusion cloning into p3XFLAG-CMV-10 vector

---

## Key Results

### Pathway Analysis
**Upregulated:**
- MYC targets
- DNA repair
- Cell cycle (G2M checkpoint, E2F targets)

**Downregulated:**
- Metabolic pathways (e.g., adipogenesis)
- Interferon response (later stages)

### Transposable Element Activation
- Strong upregulation of ERV1 family
- Moderate changes in ERVL and LINE-1
- HERV-int (ERV1) most upregulated

### Intragenic TE Enrichment
- Many upregulated TEs located within gene bodies
- Suggests roles in gene regulation and splicing

### ChIP-seq Findings
- ~154 PRMT5 binding peaks identified
- Binding observed on LINE-1 elements (especially older families)

### Locus-Level Insights
- Most loci unchanged
- Subset significantly altered
- Example: L1FLnI_10p11.22u overlapping ITGB1

### Interferon Response
- No consistent activation observed
- Some IFN-related genes detected but not supported by GSEA

---

## Conclusion
PRMT5 exhibits a dual regulatory role in epigenetic silencing and immune modulation. Its inhibition leads to selective activation of transposable elements and significant changes in cell-cycle and DNA repair pathways, without consistent interferon activation.
