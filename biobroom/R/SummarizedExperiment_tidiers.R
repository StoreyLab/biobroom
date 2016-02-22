#' Tidying methods for Biobase's SummarizedExperiment objects
#'
#' @param x SummarizedExperiment object
#' @param addPheno whether columns should be included in the tidied output
#' for those in the SummarizedExperiment colData
#' @param assay Which assay to return as the \code{value} column. Defaults to
#'   \code{assays(x)[[1L]]}
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
#'     data(airway)
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
tidy.RangedSummarizedExperiment <- function(x, addPheno=FALSE,
                                            assay=SummarizedExperiment::assayNames(x)[1L],
                                            ...) {
    if (!assay %in% SummarizedExperiment::assayNames(x)) {
        stop("Invalid assay specified: ", assay)
    }
    expressions <- SummarizedExperiment::assays(x)[[assay]] %>%
        fix_data_frame(newcol="gene")
    ret <- expressions %>%
        tidyr::gather(sample, value, -gene) %>%
        dplyr::mutate(sample=as.character(sample))

    if (addPheno) {
        pdat <- as.data.frame(colData(x))
        rownames(pdat) <- colnames(x)
        ret <- cbind(
            ret[, c('gene', 'sample')],
            pdat[ret$sample,,drop=FALSE],
            value=ret$value) %>%
            unrowname
    }
    finish(ret)
}

### Maybe add augment to add regions of interest data
