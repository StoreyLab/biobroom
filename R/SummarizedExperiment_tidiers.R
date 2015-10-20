#' Tidying methods for Biobase's SummarizedExperiment objects
#'
#' @param x SummarizedExperiment object
#' @param addPheno whether columns should be included in the tidied output
#' for those in the SummarizedExperiment colData
#' @param ... extra arguments (not used)
#'
#' @details \code{addPheno=TRUE} adds columns that are redundant (since they
#' add per-sample information to a per-sample-per-gene data frame), but that
#' are useful for some kinds of graphs and analyses.
#'
#' @name SummarizedExperiment_tidiers
#'
#' @examples
#' if (require("SummarizedExperiment", "airway")) {
#'     library(airway)
#'
#'     se <- airway
#'     tidy(se)
#' }
#'
#' @return \code{tidy} returns a data frame with one row per gene-sample
#' combination, with columns
#'   \item{gene}{gene name}
#'   \item{sample}{sample name (from column names)}
#'   \item{value}{expressions}
#'
#' If \code{addPheno} is TRUE then information from colData
#' is added.
#'
#' @S3method tidy RangedSummarizedExperiment
#' @export tidy.RangedSummarizedExperiment
tidy.RangedSummarizedExperiment <- function(x, addPheno=FALSE, ...) {
    expressions <- fix_data_frame(SummarizedExperiment::assays(x)$counts, newcol="gene")
    ret <- expressions %>% tidyr::gather(sample.id, value, -gene)

    if (addPheno) {
        pdat <- as.data.frame(colData(x))
        ret <- unrowname(as.data.frame(cbind(gene=ret$gene,
                                             pdat[ret$sample.id, ],
                                             value=ret$value)))
    }
    finish(ret)
}

### Maybe add augment to add regions of interest data
