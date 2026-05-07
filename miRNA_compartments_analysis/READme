  # miRNA-seq Differential Expression Analysis in Pediatric Burkitt Lymphoma

## Overview

This project analyzes miRNA-seq data from pediatric Burkitt lymphoma (pBL) across three biological compartments: tumor tissue (Ti), exosomes (Exo), and cell-free plasma (CF). The study integrates differential expression analysis, biological interpretation, and machine learning to characterize compartment-specific miRNA signatures and processing states.

**Total samples:** 45 (15 treatment-naïve patients, matched across compartments, n=15 per compartment)

---

## Key Findings

| Finding | Detail |
|---------|--------|
| Distinct compartment signatures | PCA shows clear separation between Ti, Exo, and CF |
| Mature miRNA enrichment | Tumor tissue |
| Precursor/processing fragments | Enriched in exosomes |
| Systemic signal | miR-122-5p (liver) in cell-free plasma |
| Classification accuracy | Random Forest: 97.8% |
| Validation | Permutation testing confirmed biological signal |

---

## Abstract

Pediatric Burkitt lymphoma (pBL) is an aggressive malignancy where circulating miRNAs represent promising biomarkers. However, their distribution across biological compartments and their relationship to processing states remain incompletely understood.

We analyzed miRNA-seq data from 45 samples derived from 15 treatment-naïve pBL patients, including matched tumor tissue, exosomal, and cell-free plasma fractions. After quality control and filtering, 507 high-confidence miRNAs were retained from 2,982 detected miRNAs. Differential expression analysis was performed using DESeq2 to identify compartment-specific miRNA signatures.

Principal component analysis showed clear separation between compartments, indicating strong compartment-specific expression patterns. Tumor tissue was enriched for mature miRNAs, while exosomes showed higher representation of reads mapping to pre-miRNA-derived fragments. Cell-free plasma exhibited mixed systemic signals, including elevated miR-122-5p.

To evaluate whether miRNA profiles encode compartment identity, a Random Forest machine learning model was trained on normalized expression data. The model achieved 97.8% mean accuracy** in distinguishing compartments. Performance dropped significantly under label permutation, confirming that classification was driven by biologically meaningful structure.



## Methods Summary

| Step | Tool/Method | Details |
|
| Quality control | FastQC, MultiQC | Raw read quality assessment |
| Read processing | fastp | Adapter trimming, 18-40 nt filtering, Phred ≥20 |
| Alignment | Bowtie (miRDeep2 mapper.pl) | GRCh38 reference genome |
| Quantification | miRDeep2 quantifier.pl | miRBase v22.1 (mature + precursor annotations) |
| Filtering | Custom R script | ≥5 reads in ≥3 samples (2,982 → 507 miRNAs) |
| Differential expression | DESeq2 | ~ compartment design, FDR < 0.1, |log2FC| > 1 |
| Visualization | ggplot2, pheatmap, EnhancedVolcano | PCA, heatmaps, volcano plots, boxplots |
| Machine learning | scikit-learn (Python) | Random Forest classifier + permutation test |

---

## Repository Structure

├── data/
│ ├── raw/ # Raw FASTQ files (see SRA accession)
│ ├── processed/ # Count matrices and metadata
│ └── reference/ # miRBase v22.1 annotations
│
├── scripts/
│ ├── 01_qc_filtering.R # Quality control and filtering
│ ├── 02_deseq2_analysis.R # Differential expression
│ ├── 03_visualization.R # PCA, heatmaps, volcano plots
│ ├── 04_compartment_plots.R # Compartment-wise comparisons
│ └── 05_machine_learning.py # Random Forest classifier
│
├── workflow/
│ └── nextflow_pipeline.nf # Reproducible Nextflow pipeline
│
├── results/
│ ├── figures/ # PCA, volcano, heatmaps, boxplots
│ ├── tables/ # Differential expression results
│ └── models/ # Trained Random Forest model
│
├── doc/
│ └── supplementary/ # QC and alignment reports
│
├── README.md
├── requirements_R.txt
├── requirements_python.txt
└── LICENSE


---

## Requirements

### R (≥ 4.2)

```r
install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr", "ggpubr"))
BiocManager::install(c("DESeq2", "VennDiagram"))

pip install scikit-learn pandas numpy matplotlib seaborn


Results
Key Outputs
Output	Location
PCA plot	results/figures/pca_compartment.pdf
Volcano plots	results/figures/volcano_Ti_vs_Exo.pdf
Heatmap	results/figures/heatmap_differential_mirnas.pdf
DE tables	results/tables/degs_Ti_vs_Exo.csv
RF results	results/models/rf_results.json
Processing State Distribution
Compartment	Mature miRNAs	Precursor fragments
Tumor tissue	High	Low
Exosomes	Moderate	High
Cell-free plasma	Mixed	Mixed
Machine Learning Results
Metric	Value
Classifier	Random Forest
Cross-validation	Stratified 5-fold
Mean accuracy	97.8%
Permutation test	p < 0.001
Top features	miRNAs distinguishing Ti, Exo, and CF
Data Availability
Raw sequencing data: NCBI SRA PRJNA1006509

Processed count matrices: data/processed/

## Citation:
 
@misc{Muhammad Jasim 2026,
  author = {Muhammad Jasim},
  title = {Compartment-resolved miRNA profiling in pediatric Burkitt lymphoma},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/jasmu06/pbl-mirna-compartment-analysis}
}
Contact
Author: Muhammad Jasim
Background: Independent Researcher | Genetics
Email: muhammadjasim.mj@gmail.com

License
MIT License

 

 



 

 
 
 
