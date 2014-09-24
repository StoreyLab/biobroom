#' Tidying methods for Biobase's DESeqDataSet objects
#'
#' @param x ExpressionSet object
#' @param colData whether colData should be included in the tidied output
#' for those in the DESeqDataSet object
#' @param ... extra arguments (not used)
#'
#' @details \code{colDat=TRUE} adds columns that are redundant (since they
#' add per-sample information to a per-sample-per-gene data frame), but that
#' are useful for some kinds of graphs and analyses.
#'
#' @name DESeqDataSet_tidiers
#'
#' @examples
#' # Nothing yet
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
