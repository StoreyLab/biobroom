#' Tidying methods for DESeq2 DESeqDataSet objects
#'
#' @param x DESeqDataSet object
#' @param colData whether colData should be included in the tidied output
#' for those in the DESeqDataSet object
#' @param ... extra arguments (not used)
#'
#' @details \code{colDat=TRUE} adds covariates from colData to the data frame.
#'
#' @name DESeqDataSet_tidiers
#'
#' @importClassesFrom DESeq2 DESeqDataSet
#'
#' @examples
#' # From DESeq2 documentation
#'
#' library("DESeq2")
#' dds <- makeExampleDESeqDataSet()
#'
#' xx <- tidy(dds)
#'
#' # With covariates
#' xx <- tidy(dds, colData=TRUE)
#'
#' # Filter data for counts >20
#' xx %>% filter(counts > 20)
#'
#' @method tidy DESeqDataSet
#'
#' @export
setMethod("tidy", "DESeqDataSet", function(x, colData=FALSE, ...) {
    expressions <- fix_data_frame(counts(x), newcol="gene")
    ret <- expressions %>% gather(sample.id, value, -gene)

    if (colData) {
        cdat <- data.frame(colData(x))
        ret <- unrowname(as.data.frame(cbind(gene=ret$gene,
                                             cdat[ret$sample.id, ],
                                             value=ret$value)))
    }
    tbl_df(ret)
})
