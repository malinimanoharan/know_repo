---
title: Analysis of RNA Sequencing data
authors:
- Malini Manoharan
tags:
- knowledge
- example
created_at: 2020-12-08 00:00:00
updated_at: 2020-12-08 22:44:13.560320
tldr: This is short description of the content and findings of the post.
---

```r
knitr::opts_chunk$set(echo = TRUE)
#Read the count count file
tcell_data <- read.csv("data/Exp6_Exp3_Counts_Test.csv", sep=",", header= TRUE,row.names = 1)
#Read the meta data file
tcell_meta <- read.csv("data/Metadata.csv", sep=",", header= TRUE)
#Deseq2 for data normalization
library(DESeq2)
```

```
## Error in library(DESeq2): there is no package called 'DESeq2'
```

```r
library(pheatmap)
```

```
## Error in library(pheatmap): there is no package called 'pheatmap'
```

```r
library(ggplot2)
```

```
## Error in library(ggplot2): there is no package called 'ggplot2'
```

```r
library(cowplot)
```

```
## Error in library(cowplot): there is no package called 'cowplot'
```

The RNAseq Data provided is a time series data of T-cell activation across different donors in two different experiments.

## Data Normalization
All RNASeq analsis begins with a sample level QC using Principal Component Analysis (PCA) and hierarchical clustering methods. It helps us understand if the samples being analyed 
in the following aspects.

    Which samples are similar to each other, which are different?
    Does this fit to the expectation from the experiment’s design?
    What are the major sources of variation in the dataset?
    
sample-level QC allows us to see how well replicates cluster together, as well as, observe whether our experimental condition represents the major source of variation in the data. Performing sample-level QC can also identify any sample outliers, which may need to be explored further to determine whether they need to be removed prior to DE analysis.Normalized count data is used to compute the PCA. Normalization in this case is done using Deseq2. 

The SummarizedExperiment class is used to store rectangular matrices of experimental results, which are commonly produced by sequencing and microarray experiments. Each object stores observations of one or more samples, along with additional meta-data describing both the observations (features) and samples (phenotypes).



```r
dds <- DESeqDataSetFromMatrix(countData = tcell_data, colData = tcell_meta, design = ~ Time)
```

```
## Error in DESeqDataSetFromMatrix(countData = tcell_data, colData = tcell_meta, : could not find function "DESeqDataSetFromMatrix"
```

```r
dds <- estimateSizeFactors(dds)
```

```
## Error in estimateSizeFactors(dds): could not find function "estimateSizeFactors"
```

```r
normalized_counts <- counts(dds, normalized=TRUE)
```

```
## Error in counts(dds, normalized = TRUE): could not find function "counts"
```

```r
se_normalized <- SummarizedExperiment(assays=list(counts=normalized_counts), colData=tcell_meta)
```

```
## Error in SummarizedExperiment(assays = list(counts = normalized_counts), : could not find function "SummarizedExperiment"
```

```r
dt_normalized <- DESeqTransform(se_normalized)
```

```
## Error in DESeqTransform(se_normalized): could not find function "DESeqTransform"
```
## Principal Component Analysis (PCA)

Principal Component Analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset (dimensionality reduction).If two samples have similar levels of expression for the genes that contribute significantly to the variation represented by PC1, they will be plotted close together on the PC1 axis. Therefore, we would expect that replicates to have similar scores (since the same genes are changing) and cluster together on PC1 and/or PC2, and the samples from different time points/experiments to have different scores. 


```r
plotPCA(dt_normalized, intgroup = "Batch", returnData = FALSE)
```

```
## Error in plotPCA(dt_normalized, intgroup = "Batch", returnData = FALSE): could not find function "plotPCA"
```

```r
norm_dis <- dist(t(normalized_counts), method = "euclidean")
```

```
## Error in t(normalized_counts): object 'normalized_counts' not found
```

## Heirarchical Clustering


```r
norm_hc <- hclust(norm_dis, method = "ward.D")
```

```
## Error in hclust(norm_dis, method = "ward.D"): object 'norm_dis' not found
```

```r
plot(norm_hc,cex=0.5)
```

```
## Error in plot(norm_hc, cex = 0.5): object 'norm_hc' not found
```
### Observations

The PCA plot and the heirarchical clustering suggests that samples from the different experiments are clusted together and therefore is indicative of batch effects coming from the different experimental setup.  

It was also observed that the pca does not change with the top 500/5000 genes suggesting a less noisy data. This can be reconfimed with the correlation heatmap (below) which suggests that the expression of the samples across the two experiments are highly correlated with respect to the different time points. 

```r
p1<-plotPCA(dt_normalized, intgroup = "Batch", ntop = 500, returnData = FALSE)
```

```
## Error in plotPCA(dt_normalized, intgroup = "Batch", ntop = 500, returnData = FALSE): could not find function "plotPCA"
```

```r
p2<-plotPCA(dt_normalized, intgroup = "Batch", ntop = 5000, returnData = FALSE)
```

```
## Error in plotPCA(dt_normalized, intgroup = "Batch", ntop = 5000, returnData = FALSE): could not find function "plotPCA"
```

```r
plot_grid(p1, p2)
```

```
## Error in plot_grid(p1, p2): could not find function "plot_grid"
```


```r
norm_mat <- assay(dt_normalized)
```

```
## Error in assay(dt_normalized): could not find function "assay"
```

```r
norm_cor <- cor(norm_mat)
```

```
## Error in is.data.frame(x): object 'norm_mat' not found
```

```r
fontsize_row = 10 - nrow(norm_cor) / 15
```

```
## Error in nrow(norm_cor): object 'norm_cor' not found
```

```r
pheatmap(norm_cor,fontsize_row=fontsize_row,fontsize_col = fontsize_row)
```

```
## Error in pheatmap(norm_cor, fontsize_row = fontsize_row, fontsize_col = fontsize_row): could not find function "pheatmap"
```
## Batch Correction

The transformations implemented in DESeq2, vst and rlog, compute a variance stabilizing transformation which is roughly similar to putting the data on the log2 scale, while also dealing with the sampling variability of low counts. It uses the design formula to calculate the within-group variability.However it does not use the design to remove variation in the data. 

It is possible to visualize the transformed data with batch variation removed, using the removeBatchEffect function from limma. This simply removes any shifts in the log2-scale expression data that can be explained by batch. The paradigm for this operation for designs with balanced batches would be:

The above approach seemed to work better for the current dataset compared  to combat from the SVA package


```r
dds <- DESeqDataSetFromMatrix(countData=tcell_data, colData=tcell_meta, design = ~ Batch + Time)
```

```
## Error in DESeqDataSetFromMatrix(countData = tcell_data, colData = tcell_meta, : could not find function "DESeqDataSetFromMatrix"
```

```r
dds <- DESeq(dds)
```

```
## Error in DESeq(dds): could not find function "DESeq"
```

```r
vsd <- vst(dds, blind=FALSE)
```

```
## Error in vst(dds, blind = FALSE): could not find function "vst"
```

```r
mat <- assay(vsd)
```

```
## Error in assay(vsd): could not find function "assay"
```

```r
mat <- limma::removeBatchEffect(mat, vsd$Batch)
```

```
## Error in loadNamespace(name): there is no package called 'limma'
```

```r
assay(vsd) <- mat
```

```
## Error in eval(expr, envir, enclos): object 'mat' not found
```

```r
counts_batch_corrected <- assay(vsd)
```

```
## Error in assay(vsd): could not find function "assay"
```

```r
se_normalized <- SummarizedExperiment(assays=list(counts=counts_batch_corrected), colData=tcell_meta)
```

```
## Error in SummarizedExperiment(assays = list(counts = counts_batch_corrected), : could not find function "SummarizedExperiment"
```

```r
dt_normalized <- DESeqTransform(se_normalized)
```

```
## Error in DESeqTransform(se_normalized): could not find function "DESeqTransform"
```
### Observations

The plots below indicate that the batch effects have been nullified.Clustering of the samples is now based on the experimental design. And the data can be compared across different time points capturing information across multiple donors. 


```r
plotPCA(dt_normalized, intgroup = "Time", returnData = FALSE)
```

```
## Error in plotPCA(dt_normalized, intgroup = "Time", returnData = FALSE): could not find function "plotPCA"
```


```r
combat_dis <- dist(t(counts_batch_corrected), method = "euclidean")
```

```
## Error in t(counts_batch_corrected): object 'counts_batch_corrected' not found
```

```r
combat_hc <- hclust(combat_dis, method = "ward")
```

```
## Error in hclust(combat_dis, method = "ward"): object 'combat_dis' not found
```

```r
plot(combat_hc,cex=0.5)
```

```
## Error in plot(combat_hc, cex = 0.5): object 'combat_hc' not found
```

## References

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8. 

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47. doi: 10.1093/nar/gkv007.