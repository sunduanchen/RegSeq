# Introduction
RegSeq package contains the proposed RegSeq algorithm (function `Reg_Finding`), which is designed to perform the regulon analysis for RNA-Seq data, especially in small sample scenario. RegSeq can identify the key regulons that are most responsible for the difference in two groups samples. In this tutorial, we use several examples to help users executing RegSeq in real applications. 

We first load the required package.
```{r}
library("RegSeq")
options(scipen = 0)
```

# RegSeq implementation example
## Input data preparation
There are two necessary input data sources of RegSeq. The first one is a gene expression signature and the second one is a regulon network. The gene expression signature is a gene-wise vector that each gene has a score, which represents distinctive phenotypes or treatments of comparing two groups of samples. Actually, any computational method that generates a quantitative measurement of difference between the groups is applicable for obtaining this vector. Here, we use the function `get_signature` to compute the expression signature.

We load the processed esophageal carcinoma (ESCA) raw counts data from TCGA.
```{r}
load('ESCA_expression.RData')
```

In this expression matrix, each row is a gene and each column is a sample. The dimensions of ESCA dataset are
```{r}
dim(exprs)
```
which indicate there are totally 14788 genes and 16 samples. Notably, this data only keep the paired samples that diagnosed with 'normal' and 'tumor'. The detailed information of these samples are:
```{r}
colnames(exprs)
```
which show that the first 8 samples are normal samples and the latter 8 samples are matched tumor samples.

Given an above RNA-Seq raw counts data with different sample labels, we can obtain the gene signature by using the compiled function `get_signature`. We set the following function parameters:
```{r}
p <- ncol(exprs)/2
rawCounts <- exprs
groupA <- 1:p
groupB <- (p+1):(2*p)
```

Then, the signature can be got by:
```{r}
signature <- get_signature(rawCounts, groupA, groupB)
```

The scores of first 5 genes are:
```{r}
signature[1:5]
```
The positive score represents highly expressed in `groupA` and negative score represents highly expressed in `groupB`, respectively. `get_signature` directly uses the differential expression analysis results in DESeq2 package. Again, users can perform any method (e.g. fold change, Student's t-test, Mann-Whitney U test) to obtain the signature and directly use the derived signature as model input.

The second part of model input is a regulon network. `Reg_Finding` use the list object of class `regulon`, where each element represent a transcriptional regulator (transcription factor) and contains two vectors: (1) a named numeric vector indicating the mode of regulation for each target gene, whose ID is indicated by the names attribute of the vector. (2) a numeric vector indicating the confidence score for the TF-target interaction. The regulon network can be generated from networks reverse engineered with the ARACNe algorithm. We now load the regulon network built from ESCA dataset.
```{r}
load('ESCA_regulon.RData')
```

This regulon network contains the following number of transcription factors:
```{r}
length(regulon)
```


## Perform regulon analysis via RegSeq
The regulon analysis can be achieved by using the folloing code:
```{r}
set.seed(11)
result <- Reg_Finding(signature, regulon, Permutation = 10, sort = TRUE)
```

In the implementation trial, the seed can obtain a reproducible result. The algorithm can automatically print out the current process. The final output is saved in a `data.frame` table with the following elements:
```{r}
head(result)
```

In this result, each row is a transcription factor and the rows are sorted in `FDR` ascending order, due to the parameter `sort = TRUE`. The column names `Positive` and `Negative` are the number of positive and negative target genes of current transcription factor, respectively. `Target` is the total number of target genes. `Reg_Score` is the regulon score computed by RegSeq and `FDR` represents the FDR value of absolute-based GSEA. The final condition of regulon is shown in `Direction`, which is a indicator to show whether the significant regulon is being activated (`ACTIVATED`) or repressed (`REPRESSED`).

Besides, the users can set other parameters to meet different needs. For example, we want find the significant regulon (FDR cutoff < 0.05) with size range from 100 to 300 and save the final result to a file named 'result.txt'. This request can be achieved by using the following code:

```{r}
set.seed(123)
result <- Reg_Finding(signature, regulon, 100, 300, cutoff = 0.05, outfile = 'result.txt')
head(result)
```

# Reference
Duanchen Sun and Zheng Xia (2018): RegSeq: A novel GSEA-based algorithm to perform the regulon analysis on RNA-Seq data with small sample replicates.
