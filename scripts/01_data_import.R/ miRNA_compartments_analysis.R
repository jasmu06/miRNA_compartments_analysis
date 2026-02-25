# ============================================================
# miRNA-seq Differential Expression Workflow (DESeq2)
# Compartments: exo / ti / cf
# Author: Your Name
# ============================================================

# -------------------------------
# 1. Load Required Libraries
# -------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(RColorBrewer)
  library(ggrepel)
  library(cowplot)
  library(VennDiagram)
  library(ggtern)
})

# -------------------------------
# 2. Input Files
# -------------------------------
counts_file   <- "data/processed/miRNA_count_matrix.csv"
metadata_file <- "data/processed/metadata_mirna.csv"

counts <- read.csv(counts_file,
                   row.names = 1,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)

metadata <- read.csv(metadata_file,
                     row.names = 1,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

metadata <- metadata[, "Compartment", drop = FALSE]

# -------------------------------
# 3. Data Preparation
# -------------------------------
counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))
counts_mat <- as.matrix(counts)
counts_mat[is.na(counts_mat)] <- 0

metadata <- metadata[colnames(counts_mat), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = metadata,
  design    = ~ Compartment
)

dds <- dds[rowSums(counts(dds)) > 0, ]

# Use poscounts for sparse miRNA data
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, sfType = "poscounts")

# -------------------------------
# 4. Pairwise Differential Expression
# -------------------------------
res_exo_vs_cf <- results(dds, contrast = c("Compartment", "exo", "cf"))
res_exo_vs_ti <- results(dds, contrast = c("Compartment", "exo", "ti"))
res_ti_vs_cf  <- results(dds, contrast = c("Compartment", "ti",  "cf"))

save_results <- function(res, filename) {
  res <- as.data.frame(res)
  res <- res[order(res$padj), ]
  write.csv(res, filename)
}

dir.create("results", showWarnings = FALSE)

save_results(res_exo_vs_cf, "results/exo_vs_cf.csv")
save_results(res_exo_vs_ti, "results/exo_vs_ti.csv")
save_results(res_ti_vs_cf,  "results/ti_vs_cf.csv")

# -------------------------------
# 5. PCA Analysis
# -------------------------------
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vsd_mat <- assay(vsd)

pca <- prcomp(t(vsd_mat))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

pca_df <- data.frame(
  Sample = colnames(vsd),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Compartment = colData(vsd)$Compartment
)

dir.create("figures", showWarnings = FALSE)

p <- ggplot(pca_df, aes(PC1, PC2, color = Compartment)) +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Global miRNA Variation Across Compartments",
    x = paste0("PC1: ", round(percentVar[1],1), "%"),
    y = paste0("PC2: ", round(percentVar[2],1), "%")
  )

ggsave("figures/PCA_miRNA.png", p, width = 8, height = 6, dpi = 300)

# -------------------------------
# 6. Heatmap of Top Significant miRNAs
# -------------------------------
sig_list <- unique(c(
  rownames(res_exo_vs_cf)[which(res_exo_vs_cf$padj < 0.05)][1:20],
  rownames(res_exo_vs_ti)[which(res_exo_vs_ti$padj < 0.05)][1:20],
  rownames(res_ti_vs_cf)[which(res_ti_vs_cf$padj < 0.05)][1:20]
))

mat <- assay(vsd)[sig_list, ]
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

annotation_col <- data.frame(
  Compartment = colData(vsd)$Compartment
)
rownames(annotation_col) <- colnames(mat_scaled)

pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Top Differential miRNAs"
)

# -------------------------------
# 7. Mature vs Precursor Classification
# -------------------------------
norm_counts <- counts(dds, normalized = TRUE)

miRNA_type <- data.frame(
  miRNA = rownames(norm_counts),
  Type = ifelse(grepl("-5p|-3p$", rownames(norm_counts)),
                "Mature", "Precursor")
)

counts_long <- as.data.frame(norm_counts) %>%
  rownames_to_column("miRNA") %>%
  pivot_longer(-miRNA, names_to = "Sample", values_to = "Counts") %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample") %>%
  left_join(miRNA_type, by = "miRNA")

counts_summary <- counts_long %>%
  group_by(Compartment, Type) %>%
  summarise(TotalCounts = sum(Counts), .groups = "drop")

counts_summary$logCounts <- log2(counts_summary$TotalCounts + 1)

ggplot(counts_summary,
       aes(Compartment, logCounts, fill = Compartment)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Type) +
  theme_classic()

# -------------------------------
# 8. Ternary Plot (Compartment Distribution)
# -------------------------------
exo_samples <- colnames(norm_counts)[grep("_exo", colnames(norm_counts))]
ti_samples  <- colnames(norm_counts)[grep("_ti",  colnames(norm_counts))]
cf_samples  <- colnames(norm_counts)[grep("_cf",  colnames(norm_counts))]

ternary_df <- data.frame(
  miRNA = rownames(norm_counts),
  Exo = rowMeans(norm_counts[, exo_samples, drop = FALSE]),
  Ti  = rowMeans(norm_counts[, ti_samples, drop = FALSE]),
  CF  = rowMeans(norm_counts[, cf_samples, drop = FALSE])
) %>%
  filter(Exo + Ti + CF > 0) %>%
  mutate(
    Total = Exo + Ti + CF,
    Exo = Exo / Total,
    Ti  = Ti  / Total,
    CF  = CF  / Total,
    Class = case_when(
      Exo >= 0.6 ~ "Exo-enriched",
      Ti  >= 0.6 ~ "Tissue-enriched",
      CF  >= 0.6 ~ "CF-enriched",
      TRUE       ~ "Mixed"
    )
  )

ggtern(data = ternary_df,
       aes(x = Ti, y = CF, z = Exo, color = Class)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(title = "miRNA Distribution Across Compartments")
