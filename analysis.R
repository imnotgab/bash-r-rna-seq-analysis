---
title: "XYZ Project"
author: "Gabriela Rybacka"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

base_path <- "./"

# What is the project goal?

The goal of this project is to compare sequencing results obtained using MiSeq and NextSeq platforms.
The analysis aims to evaluate the consistency between the two and determine if technical differences influence the observed biological signal.

# Libraries

```{r libraries}
library(DESeq2)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(knitr)
library(ggplot2)
```

# Data loading

```{r data-loading}
data_SE <- read.table(paste0(base_path, "counts_SE.txt"),
                      header = TRUE, skip = 1)

data_PE <- read.table(paste0(base_path, "counts_PE.txt"),
                      header = TRUE, skip = 1)

data_main <- cbind(data_PE[, 7:10], data_SE[, 7:10])
rownames(data_main) <- data_PE$Geneid

colnames(data_main) <- c("SRR3191542", "SRR3191543", "SRR3191544", "SRR3191545",
  "SRR3194428", "SRR3194429", "SRR3194430", "SRR3194431")

head(data_main)
```

# Differential Expression Analysis with DESeq2

```{r deseq}
samples <- colnames(data_main)

condition <- factor(rep(c("Mock", "Zika", "Mock", "Zika"), times = 2))
instrument <- factor(rep(c("MiSeq", "NextSeq"), each = 4))

colData <- data.frame(row.names = samples,
  condition = condition,
  instrument = instrument)

dds <- DESeqDataSetFromMatrix(countData = data_main,
  colData = colData,
  design = ~ instrument)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- estimateSizeFactors(dds)
log_data <- rlog(dds)
norm_data <- as.data.frame(assay(log_data))

dds <- DESeq(dds)
res <- results(dds)

summary(res)
```

Filtering low-count genes to improve reliability and reduce background noise

# Correlation - MiSeq vs. NextSeq Comparison

```{r correlation}
Miseq_vals <- unlist(norm_data[, 1:4])
NextSeq_vals <- unlist(norm_data[, 5:8])

(cor_result <- cor.test(Miseq_vals, NextSeq_vals))
```

The Pearson correlation coefficient is 0.99 with a p-value < 2.2e-16.
This indicates extremely high consistency in gene expression levels between the two sequencing platforms.

# Visualization

## Variance Heatmap

```{r heatmap-variance}
# Selecting the 1000 most variable genes for the heatmap
countVar <- apply(norm_data, 1, var)
highVar <- order(countVar, decreasing = TRUE)[1:1000]
hmDat <- as.matrix(norm_data[highVar, ])

palette <- colorRampPalette(brewer.pal(11, "PiYG"))(100)
col.inst <- ifelse(colData$instrument == "MiSeq", "magenta", "cyan")

heatmap.2(
  hmDat,
  col = rev(palette),
  trace = "none",
  scale = "row",
  labRow = FALSE,
  main = "1000 Most Variable Genes",
  ColSideColors = col.inst
)
```

## High Expression Heatmap

```{r heatmap-high-expression}
threshold_val <- 5
filtered_data <- norm_data %>% filter_all(all_vars(. >= threshold_val))

if (nrow(filtered_data) > 0) {
  annotation_col <- data.frame(
    Platform = colData$instrument,
    Condition = colData$condition
  )
  rownames(annotation_col) <- colnames(filtered_data)

  ann_colors <- list(
    Platform = c("NextSeq" = "magenta", "MiSeq" = "cyan"),
    Condition = c("Mock" = "pink", "Zika" = "violet")
  )

  pheatmap(
    as.matrix(filtered_data),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    main = "Highly Expressed Genes",
    annotation_col = annotation_col,
    annotation_colors = ann_colors
  )
}
```

Samples cluster primarily by biological condition, indicating that the biological signal is preserved regardless of the sequencing platform used.

# Principal Component Analysis (PCA)

```{r PCA}
pcaData <- plotPCA(log_data, intgroup=c("condition", "instrument"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=instrument)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("PCA - Sample Comparison")
```

The PCA plot confirms that samples cluster according to biological condition.
Technical differences between MiSeq and NextSeq have a negligible impact on the observed biological signal, which is consistent with the Pearson correlation results.
