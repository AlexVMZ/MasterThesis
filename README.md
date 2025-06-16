# README for Single-Cell RNA-Seq Analysis Pipeline
Overview
This repository contains a comprehensive pipeline for single-cell RNA-seq data analysis, from raw data processing to advanced downstream analyses. The pipeline is structured into 7 main steps, each corresponding to a specific analysis phase.

# Pipeline Structure
1. Data Preprocessing (0_Create_SCE.R)
- Creates SingleCellExperiment objects from raw 10X Genomics data
- Performs initial quality control with empty droplet removal
- Generates barcode rank plots and quality metrics
- Output: Processed SCE objects for each sample

2. Quality Control (1_QC.R)
- Performs comprehensive quality control on single-cell data
- Filters cells based on:
  Library size
  Number of detected genes
  Mitochondrial content
  Doublet detection
- Filters low-quality genes
- Generates QC plots and summary statistics

3. Normalization (2_Normalization.R)
- Implements two normalization approaches:
  1. scran: Uses deconvolution-based size factors
  2. Seurat's SCTransform: Variance-stabilizing transformation
- Output: Normalized objects for both methods

4. Gene Ranking (3_Rank.R)
- Computes gene variability metrics (CV2, Fano factor, Gini index)
- Implements two ranking methods:
  1. Rank Products: Non-parametric method
  2. SumRank: Composite ranking based on multiple metrics
- Output: Ranked gene lists and variability statistics

5. Normalization Method Comparison (5_Norm_analyses.R)
- Compares scran and SCTransform normalization results
- Analyzes differences in gene rankings
- Performs functional enrichment on differentially ranked genes
- Generates comparative visualizations

6. Meta-Analysis (6_Meta_analisis.R)
- Identifies conserved and cell-line-specific gene signatures
- Performs intersection analysis across cell lines
- Conducts functional enrichment on shared gene sets
- Generates Venn diagrams and UpSet plots

7. Gene Set Enrichment Analysis (7_GSEA_Analysis.R)
- Performs GSEA using multiple databases:
  Gene Ontology (BP, MF, CC)
  Reactome pathways
  KEGG pathways
- Generates enrichment plots and pathway visualizations

Installation and Usage
# Prerequisites
- R (≥ 4.0)
- Bioconductor (≥ 3.12)
- Required R packages:
  install.packages(c("SingleCellExperiment", "Seurat", "scran", "scater", "DropletUtils", 
                 "dplyr", "ggplot2", "patchwork", "RankProd", "clusterProfiler",
                 "org.Hs.eg.db", "fgsea", "VennDiagram", "UpSetR", "celda"))


# Running the Pipeline
1. Update file paths in each script to match your system
2. Run scripts in numerical order (0-7)
3. For cell line-specific analyses, modify the sample_name variable

# Output
- Each analysis step generates:
- Processed data files (RDS/CSV)
- Quality control plots
- Comparative visualizations
- Statistical summaries
- Functional enrichment results

# Customization
- Key parameters to customize:
- Quality control thresholds in 1_QC.R
- Normalization parameters in 2_Normalization.R
- Ranking method combinations in 3_Rank.R
- Enrichment analysis parameters in 7_GSEA_Analysis.R

# License
This project is licensed under the MIT License.

# Contact
For questions or issues, please contact the author.
