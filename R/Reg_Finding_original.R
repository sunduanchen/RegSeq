#' Onetailed: Absolute gene permuting GSEA code.
Onetailed = function(regulon, tvalue, min, max, nPerm){
    .Call("OnetailedGSEA", regulon, tvalue, min, max, nPerm)
}

#' Twotailed: Reverse positive/negative target gene permuting GSEA.
Twotailed = function(regulon, tvalue, min, max, nPerm){
    .Call("TwotailedGSEA", regulon, tvalue, min, max, nPerm)
}

#' Reg_Finding
#'
#' Reg_Finding is a novel proposed RNA-Seq tool, which can identify the key regulons that most responsible for the difference in
#' two group samples. Reg_Finding is especially designed for RNA-Seq data with small sample replicates.
#'
#' @param rawCounts The RNA-Seq raw counts expression matrix of related dataset. Each row represents
#' a gene and each column represents a sample.
#' @param regulon Regulon object contains the information of transcription factor and target genes, which can be
#' generated from networks reverse engineered with the ARACNe algorithm.
#' @param groupA Column indexes of the first group samples. The contrast of two groups in DESeq2 is groupA vs groupB.
#' @param groupB Column indexes of the second group samples.
#' @param minsize The minimum size (number of target genes) of output regulon. Only regulons with size larger than
#' \code{minsize} are considered. The default value is 10.
#' @param maxsize The maximum size (number of target genes) of output regulon. Only regulons with size less than
#' \code{maxnsize} are considered. The default value is 1000.
#' @param Permutation The times of gene-permutation in GSEA algorithm. The default value is 1000.
#' @param cutoff The cutoff of FDR for selecting the significant TFs. The default value is 0.01.
#' @param sort Locaical variable indicate whether the final output result should be sorted in FDR ascending order.
#' @param outfile The filename of output result. The file is saved at the current working directory. If \code{outfile == NULL}, algorithm would not
#' save the final result automatically.
#'
#' @details Reg_Finding algorithm compute the gene signature using the current tools in DESeq2 package. The core idea of identifying the key regulon
#' is using the enrichment analysis on the (positive/negative) target genes of a given TF. After obtains the gene signature, two types of GSEA
#' analyses are executed. The first kind of GSEA is reverse positive/negative target genes based, which
#' can identify the important TFs according to the distribution of positive/negative genes. The second is absolute-based GSEA, which can
#' tackle the inter-gene-correlation within the target genes and reduce the false positives resulted from the gene permutation in GSEA. Finally,
#' algorithm outputs the TFs that both significant in two types of GSEA (\code{FDR <= cutoff}) and use the larger normalized enrichment score in two
#' types of GSEA as the \code{Reg_Score}.
#'
#' @return This function will return a \code{data.frame} with the following components (columns):
#'   \item{TF_Name}{The gene symbol of the output transcription factor.}
#'   \item{Positive}{The number of positive target genes of current TF.}
#'   \item{Negative}{The number of negative target genes of current TF.}
#'   \item{Target}{The total number of target genes of current TF.}
#'   \item{Reg_Score}{The regulon score computed by Reg_Finding. See details for more description.}
#'   \item{FDR}{The FDR value of absolute-based GSEA.}
#'   \item{Direction}{The indicator to show whether the significant TF is \code{ACTIVATED} or \code{REPRESSED}.}
#'
#' @references
#' Duanchen Sun and Zheng Xia (2018): RegSeq: A novel GSEA-based algorithm to perform the regulon analysis on RNA-Seq
#' data with small sample replicates.
#'
#' @examples
#' # Identifying the key transcription factors using RegSeq package.
#' data(example)
#' load(regulon_info)
#' groupA <- which(dataset$label == 'Tumor')
#' groupB <- which(dataset$label == 'Normal')
#' system.time(result <- Reg_Finding(dataset$counts, regulon_info, groupA, groupB))
#'
#' @import DESeq2
#' @import Rcpp
#' @useDynLib RegSeq
#'
#' @export
Reg_Finding_original <- function(rawCounts, regulon, groupA, groupB, minsize=10, maxsize=1000, Permutation=1000, cutoff=0.01, sort=TRUE, outfile='RegSeq_Outputs.txt'){

    # Using DESeq2 to compute the gene signature score from raw counts.
    print('STEP1: Using DESeq2 to compute the gene signature score from raw counts...')
    rawCounts  <- round(rawCounts)
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
    score <- results(dds, contrast = c('condition', 'groupA', 'groupB'))[,4]
    score[is.na(score)] <- 0
    names(score) <- genenames

    # Regulon pre-processing.
    print('STEP2: Regulon preprocessing start...')
    p <- length(regulon)
    target <- list()
    for (i in 1:p){
        index <- is.element(names(regulon[[i]]$tfmode), genenames)
        regulon[[i]]$tfmode <- regulon[[i]]$tfmode[index]
        target[[i]] <- names(regulon[[i]]$tfmode)
    }

    # Absolute gene permuting GSEA.
    print('STEP3: Run absolute gene permuting GSEA...')
    genescore_abs <- abs(score)
    genescore_abs <- sort(genescore_abs, decreasing = TRUE)
    Onetailed_result <- Onetailed(regulon, genescore_abs, min = minsize, max = maxsize, nPerm = 1000)
    rownames(Onetailed_result) <- Onetailed_result[,1]
    Absolute_table <- Onetailed_result[,-1]

    # Result outputs arrangement.
    index <- which(Absolute_table[,'FDR.Q.value'] < cutoff)
    infos <- matrix(NA, length(index), 6)
    print(length(index))
    rownames(infos) <- as.character(Onetailed_result[index,1])
    colnames(infos) <- c('Positive', 'Negative', 'Target', 'Reg_Score', 'FDR', 'Direction')
    genescore_pos <- genescore_neg <- list()
    ind <- NULL
    for (i in 1:nrow(infos)){
        ind[i] <- which(names(regulon) == rownames(infos)[i])
        positive <- names(which(regulon[[ind[i]]]$tfmode > 0))
        negative <- names(which(regulon[[ind[i]]]$tfmode < 0))
        infos[i,'Positive'] <- length(positive)
        infos[i,'Negative'] <- length(negative)
        infos[i,'Target'] <- length(positive) + length(negative)
        score_pos <- score_neg <- score
        score_pos[positive] <- (-1)*(score_pos[positive])
        score_neg[negative] <- (-1)*(score_neg[negative])
        genescore_pos[[i]] <- sort(score_pos, decreasing = TRUE)
        genescore_neg[[i]] <- sort(score_neg, decreasing = TRUE)
    }
    names(genescore_pos) <- names(genescore_neg) <- names(regulon)[ind]
    Result_table <- data.frame(infos, stringsAsFactors = FALSE)
    Result_table$FDR <- Absolute_table$FDR.Q.value[index]

    # Reverse positive/negative target gene permuting GSEA.
    print('STEP4: Run reverse positive/negative target gene permuting GSEA....')
    TwoTailed_pos <- Twotailed(regulon[ind], genescore_pos, min = minsize, max = maxsize, nPerm = Permutation)
    TwoTailed_neg <- Twotailed(regulon[ind], genescore_neg, min = minsize, max = maxsize, nPerm = Permutation)

    reg_score <- pmax(abs(TwoTailed_pos$NES), abs(TwoTailed_neg$NES), abs(Absolute_table$NES[index]))
    signs_ind <- which(TwoTailed_neg$Direction == 'REPRESSED')
    reg_score[signs_ind] <- reg_score[signs_ind]*(-1)
    Result_table$Reg_Score <- reg_score
    Result_table$Direction <- TwoTailed_neg$Direction
    index  <- which(TwoTailed_pos$Nominal.P.value < cutoff | TwoTailed_neg$Nominal.P.value < cutoff)
    Output <- Result_table[index,]

    # Outputs options.
    print('Algorithm finished!!!')
    if (sort == TRUE){
        Output <- Output[order(Output$FDR),]
    }
    if (!is.null(outfile) && is.character(outfile)){
        prints <- cbind(rownames(Output), Output)
        colnames(prints)[1] <- 'TF_Name'
        write.table(prints, file = outfile, quote = FALSE, row.names = FALSE, sep = '\t')
    }
    return(Output)
}
