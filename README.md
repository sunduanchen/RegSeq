# RegSeq
RegSeq package

# Introduction
RegSeq package contains the proposed RegSeq algorithm (function `Reg_Finding`), which is designed to perform the regulon analysis for RNA-Seq data, especially in small sample scenario. RegSeq can identify the key regulons that are most responsible for the difference in two groups samples. In this tutorial, we use several examples to help users executing RegSeq in real applications. 

We first load the required package.
```{r}
library("RegSeq")
options(scipen = 0)
```
