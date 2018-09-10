#' get_signature
#'
#' Compute the gene signature using the current tools in DESeq2 package.
#'
#' @param rawCounts The RNA-Seq raw counts expression matrix of related dataset. Each row represents a gene and each
#' column represents a sample.
#' @param groupA Column indexes of the first group samples. The contrast of two groups in DESeq2 is groupA vs groupB.
#' @param groupB Column indexes of the second group samples.
#' @param remove Logical variable that indicate whether the genes with missing adjusted p-values in DESeq2 result should be removed.
#'
#' @import DESeq2
#'
#' @export
get_signature <- function(rawCounts, groupA, groupB, remove = FALSE){

    rawCounts  <- round(rawCounts[,c(groupA,groupB)])
    genenames  <- rownames(rawCounts)
    conditions <- rep('groupA', ncol(rawCounts))
    conditions[groupB] <- 'groupB'
    samples <- colnames(rawCounts)
    colData <- cbind(samples, conditions)
    rownames(colData) <- samples
    colnames(colData) <- c('sample','condition')
    colData_new <- data.frame(colData)

    dds <- DESeqDataSetFromMatrix(countData = rawCounts, colData = colData_new, design = ~ condition)
    dds <- DESeq(dds, quiet = TRUE)
    DESeq2_table <- results(dds, contrast = c('condition', 'groupA', 'groupB'))
    if (remove == TRUE){
        keep <- which(!is.na(DESeq2_table[,'padj']))
    }else{
        keep <- 1:nrow(DESeq2_table)
    }
    signature <- DESeq2_table[keep,'stat']
    signature[is.na(signature)] <- 0
    names(signature) <- genenames[keep]

    return(signature)
}
