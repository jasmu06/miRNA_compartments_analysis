  #  COMPARTMENT-RESOLVED MIRNA PROFILING REVEALS PROCESSING STATE-DEPENDENT DISTRIBUTION IN PEDIATRIC BURKITT LYMPHOMA

## Overview

This project analyzes miRNA-seq data from pediatric Burkitt lymphoma (pBL) across three biological compartments: tumor tissue (Ti), exosomes (Exo), and cell-free plasma (CF). The study integrates differential expression analysis, biological interpretation, and machine learning to characterize compartment-specific miRNA signatures and processing states.

**Total samples:** 45 (15 treatment-naïve patients, matched across compartments, n=15 per compartment)

---
 ## Key Findings

- **Distinct compartment signatures**: PCA shows clear separation between Ti, Exo, and CF

- **Mature miRNA enrichment**: Tumor tissue

- **Precursor/processing fragments**: Enriched in exosomes

- **Systemic signal**: miR-122-5p (liver) in cell-free plasma

- **Classification accuracy**: Random Forest: 97.8%

- **Validation**: Permutation testing confirmed a biological signal

---

## Abstract

Pediatric Burkitt lymphoma (pBL) is an aggressive malignancy where circulating miRNAs represent promising biomarkers. However, their distribution across biological compartments and their relationship to processing states remain incompletely understood.

We analyzed miRNA-seq data from 45 samples derived from 15 treatment-naïve pBL patients, including matched tumor tissue, exosomal, and cell-free plasma fractions. After quality control and filtering, 507 high-confidence miRNAs were retained from 2,982 detected miRNAs. Differential expression analysis was performed using DESeq2 to identify compartment-specific miRNA signatures.

Principal component analysis showed clear separation between compartments, indicating strong compartment-specific expression patterns. Tumor tissue was enriched for mature miRNAs, while exosomes showed higher representation of reads mapping to pre-miRNA-derived fragments. Cell-free plasma exhibited mixed systemic signals, including elevated miR-122-5p.

To evaluate whether miRNA profiles encode compartment identity, a Random Forest machine learning model was trained on normalized expression data. The model achieved 97.8% mean accuracy** in distinguishing compartments. Performance dropped significantly under label permutation, confirming that classification was driven by biologically meaningful structure.



## Methods Summary

 ## 🧬 Methods

### 📥 Data Collection
- 45 samples from PRJNA1006509 (tissue, exosome, cell-free plasma; n=15 each)
- SRA Toolkit v3.2.1 + Nextflow v25.04.4.5950

- **miRDeep2**: GRCh38 alignment (Bowtie) + miRBase v22.1 annotation
### 🔧 Preprocessing
- **fastp v0.23.2**: Adapter trimming, size selection (18-40 nt), quality filter (Phred≥20)

### 🧮 Quantification

### 🧹 Filtering
- 2,982 → **507 miRNAs** (≥5 reads in ≥3 samples)

### 📊 Differential Expression
- **DESeq2 v1.38.3** (R v4.2.2)
- Criteria: |log2FC| > 1, adj.p < 0.1 (FDR)

### 📈 Visualization
- ggplot2, EnhancedVolcano, heatmaps, Venn diagrams

### 🤖 Machine Learning
- Random Forest (scikit-learn) with stratified cross-validation

### 🔬 miRNA Processing Analysis
- Mature (-5p/-3p) vs precursor-associated expression patterns
---

 ## Repository Structure

| Folder/File | Description |
|-------------|-------------|
| **data/** | |
| ├── raw/ | Raw FASTQ files (see SRA accession) |
| ├── processed/ | Count matrices and metadata |
| └── reference/ | miRBase v22.1 annotations |
| **scripts/** | |
| ├── 01_qc_filtering.R | Quality control & filtering |
| ├── 02_deseq2_analysis.R | Differential expression |
| ├── 03_visualization.R | PCA, heatmaps, volcano plots |
| ├── 04_compartment_plots.R | Compartment-wise comparisons |
| └── 05_machine_learning.py | Random Forest classifier |
| **workflow/** | |
| └── nextflow_pipeline.nf | Reproducible Nextflow pipeline |
| **results/** | |
| ├── figures/ | PCA, volcano, heatmaps, boxplots |
| ├── tables/ | Differential expression results |
| └── models/ | Trained Random Forest model |
| **doc/** | |
| └── supplementary/ | QC and alignment reports |
| **README.md** | Project documentation |
| **requirements_R.txt** | R package dependencies |
| **requirements_python.txt** | Python package dependencies |
| **LICENSE** | MIT License |
---


## Requirements

### R (≥ 4.2)

install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr", "ggpubr"))
BiocManager::install(c("DESeq2", "EnhancedVolcano", "VennDiagram"))

 pip install scikit-learn pandas numpy matplotlib seaborn joblib

## 📊 Results

### 📁 Key Outputs
- PCA plot: `results/figures/pca_compartment.pdf.`
- Volcano plots: `results/figures/volcano_Ti_vs_Exo.pdf.`
- Heatmap: `results/figures/heatmap_differential_mirnas.pdf`
- DE tables: `results/tables/degs_Ti_vs_Exo.csv.`
- RF results: `results/models/rf_results.json`

### 🔬 Processing State Distribution
- **Tumor tissue**: Mature (High), Precursor (Low)
- **Exosomes**: Mature (Moderate), Precursor (High)
- **Cell-free plasma**: Mature (Mixed), Precursor (Mixed)

### 🤖 Machine Learning Results
- **Classifier**: Random Forest
- **Cross-validation**: Stratified 5-fold
- **Mean accuracy**: 97.8%
- **Permutation test**: p < 0.001
- **Top features**: miRNAs distinguishing Ti, Exo, and CF

### 📦 Data Availability 
- Raw data: [SRA PRJNA1006509](https://www.ncbi.nlm.nih.gov/sra/PRJNA1006509)
- Processed matrices: `data/processed/.`
 
## Citation:
 
 @misc{muhammad2026compartment,
  author = Muhammad, Jasim,
  title = Compartment-resolved miRNA profiling in pediatric Burkitt lymphoma,
  year = 2026,
  publisher = GitHub,
  url = {https://github.com/jasmu06/pbl-mirna-compartment-analysis

Contact
Author: Muhammad Jasim
Background: Independent Researcher | Genetics
Email: muhammadjasim.mj@gmail.com

License
MIT License

 

 



 

 
 
 
