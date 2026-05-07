## ----setup, include=FALSE-----------------------------------------------------------------------------------------
# Set global options for all chunks
knitr::opts_chunk$set(
  echo = TRUE,       # show code, set FALSE if you want only output
  message = FALSE,   # hide package messages
  warning = FALSE,   # hide warnings
  cache = FALSE, 
  fig.width = 8,     # default figure width
  fig.height = 6,    # default figure height
  dpi = 150          # figure resolution
)


## -----------------------------------------------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggrepel)
library(DESeq2)
library(tidyr)
library(ggpubr)
library(tibble)
library(forcats)
library(VennDiagram)


## -----------------------------------------------------------------------------------------------------------------
counts <- read.csv("/Users/Bioinfo_data/miRNA_compartments_analysis/data/processed/miRNA_count_matrix.csv",
                   row.names = 1,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)

metadata <- read.csv("/Users/Bioinfo_data/miRNA_compartments_analysis/data/processed/metadata_mirna.csv",
                     row.names = 1,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

metadata <- metadata[, "Compartment", drop = FALSE]


## -----------------------------------------------------------------------------------------------------------------
counts[1:5, 1:5]


## -----------------------------------------------------------------------------------------------------------------
dim(counts)
head(metadata)


## -----------------------------------------------------------------------------------------------------------------
counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))
counts_mat <- as.matrix(counts)
counts_mat[is.na(counts_mat)] <- 0


## -----------------------------------------------------------------------------------------------------------------
metadata <- metadata[colnames(counts_mat), , drop = FALSE]


## -----------------------------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = metadata,
  design = ~ Compartment
)


## -----------------------------------------------------------------------------------------------------------------
keep <- rowSums(counts(dds) >= 5) >= 3

## -----------------------------------------------------------------------------------------------------------------
dds <- dds[keep, ]


## -----------------------------------------------------------------------------------------------------------------
dim(dds)


## -----------------------------------------------------------------------------------------------------------------
nrow(dds)


## -----------------------------------------------------------------------------------------------------------------
head(rownames(dds))


## -----------------------------------------------------------------------------------------------------------------
head(counts(dds)[, 1:3])


## -----------------------------------------------------------------------------------------------------------------
dds <- estimateSizeFactors(dds, type = "poscounts")


## -----------------------------------------------------------------------------------------------------------------
dds <- DESeq(dds, sfType = "poscounts")


## -----------------------------------------------------------------------------------------------------------------
levels(colData(dds)$Compartment)


## -----------------------------------------------------------------------------------------------------------------
res_exo_vs_cf <- results(dds, contrast = c("Compartment", "exo", "cf"))
res_exo_vs_ti <- results(dds, contrast = c("Compartment", "exo", "ti"))
res_ti_vs_cf  <- results(dds, contrast = c("Compartment", "ti", "cf"))


## -----------------------------------------------------------------------------------------------------------------
summary(res_exo_vs_cf)
summary(res_exo_vs_ti)
summary(res_ti_vs_cf)


## -----------------------------------------------------------------------------------------------------------------
sig_exo_vs_cf <- res_exo_vs_cf[which(res_exo_vs_cf$padj < 0.05), ]
sig_exo_vs_cf <- sig_exo_vs_cf[order(sig_exo_vs_cf$padj), ]

sig_exo_vs_ti <- res_exo_vs_ti[which(res_exo_vs_ti$padj < 0.05), ]
sig_exo_vs_ti <- sig_exo_vs_ti[order(sig_exo_vs_ti$padj), ]

sig_ti_vs_cf <- res_ti_vs_cf[which(res_ti_vs_cf$padj < 0.05), ]
sig_ti_vs_cf <- sig_ti_vs_cf[order(sig_ti_vs_cf$padj), ]


## -----------------------------------------------------------------------------------------------------------------
summary(res_exo_vs_cf)
summary(res_exo_vs_ti)
summary(res_ti_vs_cf)


## -----------------------------------------------------------------------------------------------------------------
top40_exo_vs_cf <- head(sig_exo_vs_cf, 40)
top40_exo_vs_ti <- head(sig_exo_vs_ti, 40)
top40_ti_vs_cf  <- head(sig_ti_vs_cf, 40)


## -----------------------------------------------------------------------------------------------------------------
head(sig_exo_vs_cf, 40)
head(sig_exo_vs_ti, 40)
head(sig_ti_vs_cf, 40)



## -----------------------------------------------------------------------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3


## -----------------------------------------------------------------------------------------------------------------
dds_filtered <- dds[keep, ]


## -----------------------------------------------------------------------------------------------------------------
dim(dds_filtered)

## -----------------------------------------------------------------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds_filtered, blind = FALSE)


## -----------------------------------------------------------------------------------------------------------------
vsd_mat <- assay(vsd)


## -----------------------------------------------------------------------------------------------------------------
pca <- prcomp(t(vsd_mat))


## -----------------------------------------------------------------------------------------------------------------
 percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100


## -----------------------------------------------------------------------------------------------------------------
pcaData <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Compartment = colData(vsd)$Compartment
)


## -----------------------------------------------------------------------------------------------------------------
pcaData$CompartmentLabel <- factor(
  pcaData$Compartment,
  levels = c("exo", "cf", "ti"),
  labels = c("Exo-miRNA", "CF-miRNA", "Ti-miRNA")
)

## -----------------------------------------------------------------------------------------------------------------
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = CompartmentLabel)) +
  geom_point(size = 5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Global miRNA Variation Across Compartments",
    x = paste0("PC1: ", round(percentVar[1], 1), "% variance"),
    y = paste0("PC2: ", round(percentVar[2], 1), "% variance"),
    color = "Compartment"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

print(p)
ggsave("PCA_plot.png", plot = p, width = 6, height = 5, dpi = 300)
#### Step 8:  Plotting PCA


## -----------------------------------------------------------------------------------------------------------------
top_miRNAs <- unique(c(
  rownames(sig_exo_vs_cf)[1:20],
  rownames(sig_exo_vs_ti)[1:20],
  rownames(sig_ti_vs_cf)[1:20]
))


## -----------------------------------------------------------------------------------------------------------------
mat <- assay(vsd)[top_miRNAs, ]
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0  


## -----------------------------------------------------------------------------------------------------------------
annotation_col <- data.frame(
  Compartment = colData(dds)$Compartment
)


## -----------------------------------------------------------------------------------------------------------------
rownames(annotation_col) <- colnames(mat_scaled)
annotation_col$Compartment <- factor(
  annotation_col$Compartment,
  levels = c("exo", "cf", "ti"),
  labels = c("Exo-miRNA", "CF-miRNA", "Ti-miRNA")
)


## -----------------------------------------------------------------------------------------------------------------
ann_colors <- list(
  Compartment = c(
    "Exo-miRNA" = "#1f77b4",  # blue
    "CF-miRNA"  = "#ff7f0e",  # orange
    "Ti-miRNA"  = "#2ca02c"   # green
  )
)


## -----------------------------------------------------------------------------------------------------------------
pdf("heatmap_miRNA.p", width = 7, height = 8)
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  border_color = NA,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "  Top Differentially Expressed miRNAs Across Compartments  "
)
dev.off()


## -----------------------------------------------------------------------------------------------------------------
res_exo_vs_cf <- results(dds, contrast = c("Compartment", "exo", "cf"))
res_exo_vs_ti <- results(dds, contrast = c("Compartment", "exo", "ti"))
res_cf_vs_ti <- results(dds, contrast = c("Compartment", "cf", "ti"))


## -----------------------------------------------------------------------------------------------------------------
res_exo_vs_cf_df <- as.data.frame(res_exo_vs_cf)
res_exo_vs_cf_df$miRNA <- rownames(res_exo_vs_cf_df)

## -----------------------------------------------------------------------------------------------------------------
res_exo_vs_ti_df <- as.data.frame(res_exo_vs_ti)
res_exo_vs_ti_df$miRNA <- rownames(res_exo_vs_ti_df)

## -----------------------------------------------------------------------------------------------------------------
res_cf_vs_ti_df <- as.data.frame(res_cf_vs_ti)
res_cf_vs_ti_df$miRNA <- rownames(res_cf_vs_ti_df)


## -----------------------------------------------------------------------------------------------------------------
# exo vs cf
res_exo_vs_cf_df$padj[is.na(res_exo_vs_cf_df$padj)] <- 1
res_exo_vs_cf_df$padj[res_exo_vs_cf_df$padj == 0] <- 1e-300

# exo vs ti
res_exo_vs_ti_df$padj[is.na(res_exo_vs_ti_df$padj)] <- 1
res_exo_vs_ti_df$padj[res_exo_vs_ti_df$padj == 0] <- 1e-300

# cf vs ti
res_cf_vs_ti_df$padj[is.na(res_cf_vs_ti_df$padj)] <- 1
res_cf_vs_ti_df$padj[res_cf_vs_ti_df$padj == 0] <- 1e-300


## -----------------------------------------------------------------------------------------------------------------
# exo vs cf
res_exo_vs_cf_df$Significance <- "NS"
res_exo_vs_cf_df$Significance[res_exo_vs_cf_df$padj < 0.05 & res_exo_vs_cf_df$log2FoldChange > 0] <- "Up"
res_exo_vs_cf_df$Significance[res_exo_vs_cf_df$padj < 0.05 & res_exo_vs_cf_df$log2FoldChange < 0] <- "Down"

# exo vs ti
res_exo_vs_ti_df$Significance <- "NS"
res_exo_vs_ti_df$Significance[res_exo_vs_ti_df$padj < 0.05 & res_exo_vs_ti_df$log2FoldChange > 0] <- "Up"
res_exo_vs_ti_df$Significance[res_exo_vs_ti_df$padj < 0.05 & res_exo_vs_ti_df$log2FoldChange < 0] <- "Down"

# cf vs ti
res_cf_vs_ti_df$Significance <- "NS"
res_cf_vs_ti_df$Significance[res_cf_vs_ti_df$padj < 0.05 & res_cf_vs_ti_df$log2FoldChange > 0] <- "Up"
res_cf_vs_ti_df$Significance[res_cf_vs_ti_df$padj < 0.05 & res_cf_vs_ti_df$log2FoldChange < 0] <- "Down"


## -----------------------------------------------------------------------------------------------------------------
plot_volcano <- function(res_df, title, filename = NULL) {

  top_labels <- res_df %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(10)

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    geom_text_repel(
      data = top_labels,
      aes(label = miRNA),
      size = 3.5,
      color = "black",
      max.overlaps = 20,
      box.padding = 0.4
    ) +
    theme_minimal() +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

  if (!is.null(filename)) {
    ggsave(filename, plot = p, width = 8, height = 6)
  }

  return(p)
}
#### Step 5: Plot volcano plots for each comparison

## -----------------------------------------------------------------------------------------------------------------
 if(!dir.exists("figures")) dir.create("figures", recursive = TRUE)
plot_volcano(res_exo_vs_cf_df, "Exo vs CF-miRNA Comparison", "figures/volcano_exo_vs_cf.pdf")


## -----------------------------------------------------------------------------------------------------------------
plot_volcano(res_exo_vs_ti_df, "Exo vs Ti-miRNA Comparison", "figures/volcano_exo_vs_ti.pdf")


## -----------------------------------------------------------------------------------------------------------------
plot_volcano(res_cf_vs_ti_df, "CF vs Ti-miRNA Comparison", "figures/volcano_cf_vs_ti.pdf")


## -----------------------------------------------------------------------------------------------------------------
res_exo_vs_cf_df$Comparison <- "Exo vs CF"
res_exo_vs_ti_df$Comparison <- "Exo vs Ti"
res_cf_vs_ti_df$Comparison <- "CF vs Ti"


## -----------------------------------------------------------------------------------------------------------------
combined_volcano <- rbind(
  res_exo_vs_cf_df,
  res_exo_vs_ti_df,
  res_cf_vs_ti_df
)

combined_volcano <- combined_volcano[!is.na(combined_volcano$padj), ]

combined_volcano$Comparison <- factor(
  combined_volcano$Comparison,
  levels = c("Exo vs Ti", "CF vs Ti", "Exo vs CF")
)


## -----------------------------------------------------------------------------------------------------------------
combined_volcano <- combined_volcano %>%
  mutate(
    significance = ifelse(padj < 0.05 & abs(log2FoldChange) > 1,
                          "Significant",
                          "Not Significant"),
    direction = ifelse(log2FoldChange > 0, "Up", "Down")
  )

top_hits <- combined_volcano %>%
  filter(padj < 0.05) %>%
  group_by(Comparison, direction) %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  ungroup()


## -----------------------------------------------------------------------------------------------------------------
 p <- ggplot(combined_volcano,
            aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6) +
  geom_point(data = top_hits, color = "blue", size = 2) +
  geom_text_repel(
    data = top_hits,
    aes(label = miRNA),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_color_manual(values = c(
    "Not Significant" = "grey70",
    "Significant" = "red"
  )) +
  theme_minimal() +
  labs(
    title = "Significant and Non-Significant miRNAs Across Compartments",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12)
  )

p

## -----------------------------------------------------------------------------------------------------------------


## -----------------------------------------------------------------------------------------------------------------
head(rownames(dds), 20)


## -----------------------------------------------------------------------------------------------------------------
miRNA_df <- data.frame(miRNA = rownames(dds))


## -----------------------------------------------------------------------------------------------------------------
miRNA_df$Type <- ifelse(
  grepl("-[35]p$", miRNA_df$miRNA),
  "Mature",
  "Precursor"
)


## -----------------------------------------------------------------------------------------------------------------
head(miRNA_df, 10)


## -----------------------------------------------------------------------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)


## -----------------------------------------------------------------------------------------------------------------
norm_df <- as.data.frame(norm_counts)


## -----------------------------------------------------------------------------------------------------------------
norm_df$miRNA <- rownames(norm_df)


## -----------------------------------------------------------------------------------------------------------------
head(norm_df, 10) 


## -----------------------------------------------------------------------------------------------------------------
miRNA_type <- data.frame(
  miRNA = rownames(norm_counts),
  Type = miRNA_df$Type
)


## -----------------------------------------------------------------------------------------------------------------
nrow(norm_counts)
head(miRNA_type)


## -----------------------------------------------------------------------------------------------------------------
exo_samples <- colnames(norm_counts)[grep("_exo", colnames(norm_counts))]
ti_samples  <- colnames(norm_counts)[grep("_ti", colnames(norm_counts))]
cf_samples  <- colnames(norm_counts)[grep("_cf", colnames(norm_counts))]


## -----------------------------------------------------------------------------------------------------------------
head(exo_samples)
head(ti_samples)
head(cf_samples)


## -----------------------------------------------------------------------------------------------------------------
# Function to calculate total counts per miRNA per compartment for any two groups
prepare_bar_data <- function(norm_counts, miRNA_type, group1, group2, name1, name2) {
  
  # Subset counts for selected groups
  counts_subset <- norm_counts[, c(group1, group2)]
  
  # Sum counts per miRNA
  total_counts <- data.frame(
    miRNA = rownames(counts_subset),
    Comp1 = rowSums(counts_subset[, group1, drop = FALSE]),
    Comp2 = rowSums(counts_subset[, group2, drop = FALSE])
  )
  
  # Rename columns to meaningful compartment names
  colnames(total_counts) <- c("miRNA", name1, name2)
  
  # Merge with miRNA type
  total_counts <- merge(total_counts, miRNA_type, by = "miRNA")
  
 # Pivot longer for ggplot
  counts_long <- total_counts %>%
    pivot_longer(
      cols = c(name1, name2),
      names_to = "Compartment",
      values_to = "TotalCounts"
    )
  
  return(counts_long)
}

  


## -----------------------------------------------------------------------------------------------------------------
# Exo vs Ti
data_exo_ti <- prepare_bar_data(norm_counts, miRNA_type, exo_samples, ti_samples, "Exo", "Ti")

# Exo vs CF
data_exo_cf <- prepare_bar_data(norm_counts, miRNA_type, exo_samples, cf_samples, "Exo", "CF")

# Ti vs CF
data_ti_cf  <- prepare_bar_data(norm_counts, miRNA_type, ti_samples, cf_samples, "Ti", "CF")


## -----------------------------------------------------------------------------------------------------------------
head(data_exo_ti)


## -----------------------------------------------------------------------------------------------------------------
data_exo_ti  <- data_exo_ti  %>% mutate(Comparison = "Exo_vs_Ti")
data_exo_cf  <- data_exo_cf  %>% mutate(Comparison = "Exo_vs_CF")
data_ti_cf   <- data_ti_cf   %>% mutate(Comparison = "Ti_vs_CF")


## -----------------------------------------------------------------------------------------------------------------
counts_all <- bind_rows(data_exo_ti, data_exo_cf, data_ti_cf)


## -----------------------------------------------------------------------------------------------------------------
counts_all <- counts_all %>%
  mutate(logCounts = log2(TotalCounts + 1))


## -----------------------------------------------------------------------------------------------------------------
  head(counts_all)


## -----------------------------------------------------------------------------------------------------------------
 p <- ggplot(counts_all,
            aes(x = Compartment,
                y = logCounts,
                fill = Compartment)) +
  
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  
  facet_wrap(~Type) +
  
  scale_fill_manual(values = c(
    "EXO" = "#b22222",
    "TI"  = "#2e8b57",
    "CF"  = "#1f4e79"
  )) +
  
  labs(
    title = "miRNA Expression Distribution per Compartment",
    x = "Compartment",
    y = "Log2 Normalized Counts"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"
  )


print(p)


## -----------------------------------------------------------------------------------------------------------------
exo_samples <- colnames(norm_counts)[grep("_exo", colnames(norm_counts))]
ti_samples  <- colnames(norm_counts)[grep("_ti", colnames(norm_counts))]
cf_samples  <- colnames(norm_counts)[grep("_cf", colnames(norm_counts))]


## -----------------------------------------------------------------------------------------------------------------
  head(exo_samples)
  head(ti_samples)
  head(cf_samples)


## -----------------------------------------------------------------------------------------------------------------
exo_miRNAs <- rownames(norm_counts)[rowSums(norm_counts[, exo_samples, drop=FALSE]) > 0]
  ti_miRNAs  <- rownames(norm_counts)[rowSums(norm_counts[, ti_samples, drop=FALSE]) > 0]
  cf_miRNAs  <- rownames(norm_counts)[rowSums(norm_counts[, cf_samples, drop=FALSE]) > 0]


## -----------------------------------------------------------------------------------------------------------------
length(exo_miRNAs)   
  length(ti_miRNAs)    
  length(cf_miRNAs)    


## -----------------------------------------------------------------------------------------------------------------
exo_miRNAs <- rownames(norm_counts)[rowSums(norm_counts[, exo_samples] > 5) >= 2]
ti_miRNAs  <- rownames(norm_counts)[rowSums(norm_counts[, ti_samples] > 5) >= 2]
cf_miRNAs  <- rownames(norm_counts)[rowSums(norm_counts[, cf_samples] > 5) >= 2]


## -----------------------------------------------------------------------------------------------------------------
venn.plot <- venn.diagram(
    x = list(
      Exo = exo_miRNAs,
      Ti  = ti_miRNAs,
      CF  = cf_miRNAs
    ),
    filename = NULL,
    fill = c("#b22222", "#2e8b57", "#1f4e79"),
    alpha = 0.5,
    cex = 1.2,
    cat.cex = 1.3,
    cat.pos = c(-20, 20, 0),   # keep simple
    cat.dist = c(0.05, 0.05, 0.07),  # small push for CF
    lwd = 2,
    main = "Overlap of Detected miRNAs Across Compartments",
    main.cex = 1.5
  )
   grid.draw(venn.plot)


## -----------------------------------------------------------------------------------------------------------------
# Create figures directory (check carefully: 'if' not 'f')
if(!dir.exists("figures")) dir.create("figures")

 # Create plot object
venn.plot <- venn.diagram(
  x = list(
    Exo = exo_miRNAs,
    Ti  = ti_miRNAs,
    CF  = cf_miRNAs
  ),
  filename = NULL,
  fill = c("#b22222", "#2e8b57", "#1f4e79"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.3,
  lwd = 2
)

# Save manually
pdf("figures/venn_miRNA.pdf", width = 6, height = 6)
grid.draw(venn.plot)

