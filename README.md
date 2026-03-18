# PRMT5 and Retrotransposon Regulation

**Investigating the Role of PRMT5 in Regulating Retrotransposon Expression and Interferon Response**  
SRFP 2025 Project – Indian Academy of Sciences  
Author: Vaibhav D  

---

## Overview
This project investigates how PRMT5 inhibition affects **transposable element (TE) expression** and associated **cellular pathways**, integrating **GSEA (pathway-level analysis)** with **TE-focused analysis (TEtranscripts and Telescope)**.

---

## Background
Transposable elements constitute nearly half of the human genome and are typically epigenetically silenced. PRMT5 plays a critical role in maintaining this repression. Disruption of PRMT5 can lead to:
- Reactivation of retrotransposons (ERVs, LINE-1)
- Accumulation of nucleic acids
- Potential activation of immune pathways  

This study connects **TE activation with pathway-level changes** using GSEA.

---

## Objectives
- Quantify differential expression of transposable elements upon PRMT5 inhibition  
- Perform pathway-level analysis using GSEA  
- Compare TE activation across treatments (GSK591, LLY283)  
- Distinguish TE expression at:
  - Family level (TEtranscripts)
  - Locus level (Telescope)  

---

## Methodology

### Gene Set Enrichment Analysis (GSEA)
- Dataset: GEO (GSE236930)  
- Tool: **fgsea (R)**  
- Gene sets:
  - Hallmark  
  - C2 (KEGG/Reactome)  
  - C5 (GO)  

**Approach:**
- Differential expression using DESeq2  
- Ranked gene list used for enrichment analysis  
- Focus on:
  - Interferon pathways  
  - Cell cycle regulation  
  - DNA repair pathways  

---

### Transposable Element Analysis (TEtranscripts)
- Tool: **TEtranscripts**
- Handles multi-mapped reads using EM algorithm  
- Outputs:
  - TE family-level expression  
  - Differential TE activation  

---

### Intragenic vs Intergenic TE Classification
- Tool: **BEDTools**
- Classifies TE loci as intragenic or intergenic  
- Helps identify potential regulatory interactions  

---

### Locus-Level TE Quantification (Telescope)
- Tool: **Telescope**
- Resolves **specific TE insertion sites**  
- Uses EM algorithm for locus-level assignment  

---

### RNA-seq Integration
- Alignment: STAR  
- QC: FastQC, Trim Galore  
- Differential expression: DESeq2  

---

## Key Results

### GSEA Findings
- Upregulated pathways:
  - MYC targets  
  - DNA repair  
  - Cell cycle (G2M checkpoint, E2F targets)  

- Downregulated pathways:
  - Metabolic pathways (e.g., adipogenesis)  
  - Interferon response (later timepoints)  

---

### TE Family-Level Activation
- Strong upregulation of **ERV1 family elements**
- HERV-int identified as most upregulated  
- Moderate changes in LINE-1 and ERVL  

---

### TE Genomic Context
- Enrichment of **intragenic TEs**  
- Suggests roles in gene regulation and transcription  

---

### Locus-Level Insights (Telescope)
- Majority of LINE-1 loci unchanged  
- Subset significantly altered  

**Key example:**
- L1FLnI_10p11.22u (overlapping *ITGB1*) upregulated  

---

### Comparative Analysis (LLY283 vs GSK591)
- Distinct TE activation patterns  
- Context-dependent regulation of LINE-1 loci  

---

### Interferon Response
- No consistent global activation observed  
- IFN-related genes not supported by GSEA  

---

## Conclusion
PRMT5 inhibition results in coordinated changes at both **pathway level (GSEA)** and **transposable element level**, with strong activation of ERV1 families and selective regulation of LINE-1 loci. The lack of consistent interferon activation suggests a **complex and context-dependent relationship between TE activation and immune signaling**.
