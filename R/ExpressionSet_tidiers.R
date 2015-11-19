#' Tidying methods for Biobase's ExpressionSet objects
#'
#' @param x ExpressionSet object
#' @param addPheno whether columns should be included in the tidied output
#' for those in the ExpressionSet's phenoData
#' @param assay The name of the \code{\link[Biobase]{assayDataElement}} to use
#'   as the values to tidy. Defaults to \code{assayDataElementNames(x)[1L]},
#'   which is usually equivalent to \code{exprs(x)}.
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
#' @importFrom Biobase assayDataElement assayDataElementNames
tidy.ExpressionSet <- function(x, addPheno=FALSE,
                               assay=Biobase::assayDataElementNames(x)[1L],
                               ...) {
    if (!assay %in% Biobase::assayDataElementNames(x)) {
        stop("Invalid assayDataElementName: ", assay)
    }
    expressions <- Biobase::assayDataElement(x, assay) %>%
        fix_data_frame(newcol="gene")
    ret <- expressions %>%
        tidyr::gather(sample, value, -gene) %>%
        dplyr::mutate(sample=as.character(sample))

    if (addPheno) {
        pdat <- pData(x)
        rownames(pdat) <- colnames(x)
        ret <- cbind(
            ret[, c('gene', 'sample')],
            pdat[ret$sample,,drop=FALSE],
            value=ret$value) %>%
            unrowname
    }
    finish(ret)
}
