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
#' @examples
#' library(Biobase)
#' # import ExpressionSet object
#' data(hammer)
#'
#' # Use tidy to extract genes, sample ids and measured value
#' tidy(hammer)
#' # add phenoType data
#' tidy(hammer, addPheno=TRUE)
#'
#' @return \code{tidy} returns a data frame with one row per gene-sample
#' combination, with columns
#'   \item{gene}{gene name}
#'   \item{sample}{sample name (from column names)}
#'   \item{value}{expressions on log2 scale}
#'
#' @S3method tidy ExpressionSet
#' @export tidy.ExpressionSet
tidy.ExpressionSet <- function(x, addPheno=FALSE, ...) {
    expressions <- fix_data_frame(Biobase::exprs(x), newcol="gene")
    ret <- expressions %>% tidyr::gather(sample.id, value, -gene)

    if (addPheno) {
        pdat <- pData(x)
        ret <- unrowname(as.data.frame(cbind(gene=ret$gene,
                                             pdat[as.character(ret$sample), , drop=FALSE],
                                             value=ret$value)))
    }
    finish(ret)
}
