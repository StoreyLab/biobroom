#' Tidiers for edgeR's differential expression objects
#'
#' Tidy, augment and glance methods for turning edgeR objects into tidy data
#' frames, where each row represents one observation and each column represents
#' one column.
#'
#' @name edgeR_tidiers
#'
#' @param x DGEExact object
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
#' @param alpha Confidence level to test for significance
#' @param p.adjust.method Method for adjusting p-values to determine
#' significance; can be any in p.adjust.methods
#' @param ... extra arguments (not used)
#'
#' @return \code{glance} returns one row with the columns
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
