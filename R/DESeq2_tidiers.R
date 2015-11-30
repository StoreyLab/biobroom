#' Tidying methods for DESeq2 DESeqDataSet objects
#'
#' This reshapes a DESeq2 expressionset object into a tidy format. If the
#' dataset contains hypothesis test results (p-values and estimates), this
#' summarizes one row per gene per possible contrast.
#'
#' @param x DESeqDataSet object
#' @param colData whether colData should be included in the tidied output
#' for those in the DESeqDataSet object. If dataset includes hypothesis test
#' results, this is ignored
#' @param intercept whether to include hypothesis test results from the
#' (Intercept) term. If dataset does not include hypothesis testing,
#' this is ignored
#' @param ... extra arguments (not used)
#'
#' @details \code{colDat=TRUE} adds covariates from colData to the data frame.
#'
#' @return If the dataset contains results (p-values and log2 fold changes),
#' the result is a data frame with the columns
#'   \item{term}{The contrast being tested, as given to
#'   \code{results}}
#'   \item{gene}{gene ID}
#'   \item{baseMean}{mean abundance level}
#'   \item{estimate}{estimated log2 fold change}
#'   \item{stderror}{standard error in log2 fold change estimate}
#'   \item{statistic}{test statistic}
#'   \item{p.value}{p-value}
#'   \item{p.adjusted}{adjusted p-value}
#'
#' If the dataset does not contain results (\code{DESeq} has
#' not been run on it), \code{tidy} defaults to tidying the counts in
#' the dataset:
#'   \item{gene}{gene ID}
#'   \item{sample}{sample ID}
#'   \item{count}{number of reads in this gene in this sample}
#'
#' If \code{colData = TRUE}, it also merges this with the columns present
#' in \code{colData(x)}.
#'
#' @name DESeq2_tidiers
#'
#' @examples
#'
#' # From DESeq2 documentation
#'
#' if (require("DESeq2")) {
#'     dds <- makeExampleDESeqDataSet(betaSD = 1)
#'
#'     tidy(dds)
#'     # With design included
#'     tidy(dds, colData=TRUE)
#'
#'     # add a noise confounding effect
#'     colData(dds)$noise <- rnorm(nrow(colData(dds)))
#'     design(dds) <- (~ condition + noise)
#'
#'     # perform differential expression tests
#'     ddsres <- DESeq(dds, test = "Wald")
#'     # now results are per-gene, per-term
#'     tidied <- tidy(ddsres)
#'     tidied
#'
#'     if (require("ggplot2")) {
#'         ggplot(tidied, aes(p.value)) + geom_histogram(binwidth = .05) +
#'             facet_wrap(~ term, scale = "free_y")
#'     }
#' }
#'
#' @S3method tidy DESeqDataSet
#' @export tidy.DESeqDataSet
tidy.DESeqDataSet <- function(x, colData = FALSE, intercept = FALSE, ...) {
    # try to extract the per-gene, per-coefficient information
    resnames <- DESeq2::resultsNames(x)
    if (length(resnames) > 0) {
        ret <- data.frame(term = resnames) %>%
            dplyr::group_by(term) %>%
            dplyr::do(tidy(DESeq2::results(x, name = as.character(.$term))))
        if (!intercept) {
            ret <- ret %>% dplyr::ungroup() %>% dplyr::filter(term != "(Intercept)")
        }
        return(finish(ret))
    }

    # otherwise, tidy the expression data within it
    expressions <- fix_data_frame(counts(x), newcol="gene")
    ret <- expressions %>%
        tidyr::gather(sample, count, -gene) %>%
        dplyr::mutate(sample=as.character(sample))

    if (colData) {
        cdat <- data.frame(SummarizedExperiment::colData(x), stringsAsFactors=FALSE)
        rownames(cdat) <- colnames(x)
        ret <- cbind(
            ret[, c('gene', 'sample')],
            cdat[ret$sample,,drop=FALSE],
            count=ret$count) %>%
            unrowname
    }
    finish(ret)
}


#' @rdname DESeq2_tidiers
#'
#' @S3method tidy DESeqResults
#' @export tidy.DESeqResults
tidy.DESeqResults <- function(x, ...) {
    nn <- c("baseMean", "estimate", "stderror", "statistic", "p.value", "p.adjusted")
    finish(fix_data_frame(x, newnames = nn, newcol = "gene"))
}
