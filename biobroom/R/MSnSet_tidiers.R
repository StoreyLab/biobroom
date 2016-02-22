#' Tidying methods for Biobase's ExpressionSet objects
#'
#' @param x MSnSet object
#' @param addPheno whether columns should be included in the tidied output
#' for those in the MSnSet's phenoData
#' @param ... extra arguments (not used)
#'
#' @details \code{addPheno=TRUE} adds columns that are redundant (since they
#' add per-sample information to a per-sample-per-gene data frame), but that
#' are useful for some kinds of graphs and analyses.
#'
#' @name MSnSet_tidiers
#'
#' @examples
#' if (require("MSnbase")) {
#'   library(MSnbase)
#'   # import MSnSet object
#'   data(msnset)
#'
#'   # Use tidy to extract genes, sample ids and measured value
#'   tidy(msnset)
#'   # add phenoType data
#'   tidy(msnset, addPheno=TRUE)
#' }
#' @return \code{tidy} returns a data frame with one row per gene-sample
#' combination, with columns
#'   \item{protein}{protein name}
#'   \item{sample}{sample name (from column names)}
#'   \item{value}{protein quantitation data}
#'
#' @S3method tidy MSnSet
#' @export tidy.MSnSet
tidy.MSnSet <- function(x, addPheno=FALSE, ...) {
    expressions <- fix_data_frame(Biobase::exprs(x), newcol="protein")
    ret <- expressions %>% tidyr::gather(sample.id, value, -protein)

    if (addPheno) {
        pdat <- pData(x)
        ret <- unrowname(as.data.frame(cbind(protein=ret$protein,
                                             sample = ret$sample,
                                             pdat[as.character(ret$sample), , drop=FALSE],
                                             value=ret$value)))
    }
    finish(ret)
}
