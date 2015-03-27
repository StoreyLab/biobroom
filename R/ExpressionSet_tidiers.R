#' Tidying methods for Biobase's ExpressionSet objects
#'
#' @param x ExpressionSet object
#' @param addPheno whether columns should be included in the tidied output
#' for those in the ExpressionSet's phenoData
#' @param ... extra arguments (not used)
#'
#' @details \code{addPheno=TRUE} adds columns that are redundant (since they
#' add per-sample information to a per-sample-per-gene data frame), but that
#' are useful for some kinds of graphs and analyses.
#'
#' @name ExpressionSet_tidiers
#'
#' @import broom
#' @import Biobase
#' @importFrom dplyr %>%
#'
#' @method tidy DESeqDataSet
#'
#' @examples
#'
#' # import ExpressionSet object
#' data(hammer)
#'
#' # Use tidy to extract genes, sample ids and measured value
#' tidy(hammer)
#' # add phenoType data
#' tidy(hammer, addPheno=TRUE)
#'
#' # To look at the data with counts > 30
#' library("dplyr")
#'
#' tt %>% filter(value > 30)
#' # to further filter for sample ids
#' tt %>% filter(value > 30 & sample.id =="SRX020102")
#' @export
setMethod("tidy", "ExpressionSet", function(x, addPheno=FALSE, ...) {
    expressions <- fix_data_frame(exprs(x), newcol="gene")
    ret <- expressions %>% tidyr::gather(sample.id, value, -gene)

    if (addPheno) {
        pdat <- pData(x)
        ret <- unrowname(as.data.frame(cbind(gene=ret$gene,
                                             pdat[ret$sample.id, ],
                                             value=ret$value)))
    }
    finish(ret)
})
