#' Tidiers for edgeR's differential expression objects
#'
#' Tidy, augment and glance methods for turning edgeR objects into tidy data
#' frames, where each row represents one observation and each column represents
#' one column.
#'
#' @name edgeR_tidiers
#'
#' @param x DGEExact, DGEList object
#' @param data merge data to augment. This is particularly useful when merging
#' gene names or other per-gene information. Default is NULL.
#'
#' @examples
#' if (require("edgeR")) {
#'     library(Biobase)
#'     data(hammer)
#'     hammer.counts <- exprs(hammer)[, 1:4]
#'     hammer.treatment <- phenoData(hammer)$protocol[1:4]
#'
#'     y <- DGEList(counts=hammer.counts,group=hammer.treatment)
#'     y <- calcNormFactors(y)
#'     y <- estimateCommonDisp(y)
#'     y <- estimateTagwiseDisp(y)
#'     et <- exactTest(y)
#'
#'     head(tidy(et))
#'     head(glance(et))
#' }
#'
#' @S3method tidy DGEExact
#' @export tidy.DGEExact
tidy.DGEExact <- function(x, ...) {
    ret <- fix_data_frame(x$table, c("estimate", "logCPM", "p.value"),
                          newcol = "gene")
    finish(ret)
}

#' @rdname edgeR_tidiers
#'
#' @param addSamples Merge information from samples. Default is FALSE.
#'
#' @return \code{tidy} defaults to tidying the counts in
#' the dataset:
#'   \item{gene}{gene ID}
#'   \item{sample}{sample ID}
#'   \item{count}{number of reads in this gene in this sample}
#'
#' If \code{addSamples = TRUE}, it also merges this with the sample information present
#' in \code{x$samples}.
#'
#'
#' @S3method tidy DGEList
#' @export tidy.DGEList
tidy.DGEList <- function(x, addSamples = FALSE, ...) {
    if (is.null(rownames(x$counts)) | all(rownames(x$counts) == 1:nrow(x$counts))) {
      rownames(x$counts) <- paste0("g", 1:nrow(x$counts))
    }
    expressions <- fix_data_frame(x$counts, newnames = colnames(x$counts), newcol="gene")
    ret <- expressions %>%
        tidyr::gather(sample, count, -gene) %>%
        dplyr::mutate(sample=as.character(sample))
    if (addSamples) {
        sdat <- x$samples
        rownames(sdat) <- colnames(x)
        ret <- cbind(
            ret[, c('gene', 'sample')],
            sdat[ret$sample,,drop=FALSE],
            count=ret$count) %>%
            unrowname
    }

    finish(ret)
}

#' @rdname edgeR_tidiers
#' @return \code{augment} returns per-gene information (DGEList only)
#' @S3method augment DGEList
#' @export augment.DGEList
augment.DGEList <- function(x, data = NULL, ...) {
    ret <- list()
    # other columns that can be combined (include only those that exist)
    othernames <- c("genes", "AveLogCPM", "common.dispersion", "trended.dispersion", "tagwise.dispersion")
    #newnames <- c("statistic", "p.value", "lod")
    for (i in seq_along(othernames)) {
        if (!is.null(x[[othernames[i]]]))
            ret[[othernames[i]]] <- as.numeric(x[[othernames[i]]])
    }
    if (!missing(data)) {
        ret <- cbind(as.data.frame(data), as.data.frame(ret))
    }
    if (is.null(names(ret))) stop("No columns to augment in DGEList")
    finish(as.data.frame(ret))
}

#' @rdname edgeR_tidiers
#'
#' @param alpha Confidence level to test for significance
#' @param p.adjust.method Method for adjusting p-values to determine
#' significance; can be any in p.adjust.methods
#' @param ... extra arguments (not used)
#'
#' @return \code{glance} returns one row with the columns (DGEExact only)
#'   \item{significant}{number of significant genes using desired adjustment
#'   method and confidence level}
#'   \item{comparison}{The pair of groups compared by edgeR, delimited by /}
#'
#' @S3method glance DGEExact
#' @export glance.DGEExact
glance.DGEExact <- function(x, alpha=.05, p.adjust.method="fdr", ...) {
    pvals <- x$table$PValue
    pvals.adj <- p.adjust(pvals, p.adjust.method)
    data.frame(significant=sum(pvals.adj <= alpha),
               comparison=paste(x$comparison, collapse="/"))
}
