#' Tidiers for the output of limma (linear models for microarray analysis)
#'
#' Tidy, augment, and glance methods for MArrayLM objects, which contain the
#' results of gene-wise linear models to microarray datasets. This class is
#' the output of the lmFit and eBayes functions.
#'
#' Tidying this fit computes one row per coefficient per gene, while
#' augmenting returns one row per gene, with per-gene statistics included.
#' (This is thus a rare case where the \code{augment} output has more rows
#' than the \code{tidy} output. This is a side effect of the fact that the
#' input to limma is not tidy but rather a one-row-per-gene matrix).
#'
#' @name limma_tidiers
#'
#' @param x MArrayLM object
#'
#' @return The output of tidying functions is always a data frame without
#' rownames.
#'
#' @importClassesFrom limma MArrayLM
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' if (require("limma")) {
#'     # create random data and design
#'     set.seed(2014)
#'     dat <- matrix(rnorm(1000), ncol=4)
#'     dat[, 1:2] <- dat[, 1:2] + .5  # add an effect
#'     rownames(dat) <- paste0("g", 1:nrow(dat))
#'     des <- data.frame(treatment = c("a", "a", "b", "b"),
#'                       confounding = rnorm(4))
#'
#'     lfit <- lmFit(dat, model.matrix(~ treatment + confounding, des))
#'     eb <- eBayes(lfit)
#'     head(tidy(lfit))
#'     head(tidy(eb))
#'
#'     if (require("ggplot2")) {
#'         # the tidied form puts it in an ideal form for plotting
#'         ggplot(tidy(lfit), aes(estimate)) + geom_histogram(binwidth=1) +
#'             facet_wrap(~ term)
#'         ggplot(tidy(eb), aes(p.value)) + geom_histogram(binwidth=.2) +
#'             facet_wrap(~ term)
#'     }
#' }
NULL


#' @rdname limma_tidiers
#'
#' @return \code{tidy} returns one row per gene per coefficient. It always
#' contains the columns
#'   \item{gene}{The name of the gene (extracted from the rownames of the
#'   input matrix)}
#'   \item{term}{The coefficient being estimated}
#'   \item{estimate}{The estimate of each per-gene coefficient}
#'
#' Depending on whether the object comes from \code{eBayes}, it may also
#' contain
#'   \item{statistic}{Empirical Bayes t-statistic}
#'   \item{p.value}{p-value computed from t-statistic}
#'   \item{lod}{log-of-odds score}
#'
#' @import broom
#'
#' @export
setMethod("tidy", "MArrayLM", function(x, ...) {
    coefs <- fix_data_frame(x$coefficients, newnames = colnames(x),
                            newcol="gene")
    ret <- coefs %>% tidyr::gather(term, estimate, -gene)

    # other columns that can be combined (include only those that exist)
    othernames <- c("t", "p.value", "lods")
    newnames <- c("statistic", "p.value", "lod")
    for (i in seq_along(othernames)) {
        if (!is.null(x[[othernames[i]]]))
        ret[[newnames[i]]] <- as.numeric(x[[othernames[i]]])
    }

    finish(ret)
})


#' @rdname limma_tidiers
#'
#' @param data original expression matrix; if missing, \code{augment} returns only
#' the computed per-gene statistics
#' @param ... extra arguments, not used
#'
#' @return \code{augment} returns one row per gene, containing the original
#' gene expression matrix if provided. It then adds columns containing
#' the per-gene statistics included in the MArrayLM object, each prepended
#' with a .:
#'   \item{.gene}{gene ID, obtained from the rownames of the input}
#'   \item{.sigma}{per-gene residual standard deviation}
#'   \item{.df.residual}{per-gene residual degrees of freedom}
#'
#' The following columns may also be included, depending on which have been
#' added by \code{lmFit} and \code{eBayes}:
#'   \item{.AMean}{average intensity across probes}
#'   \item{.statistic}{moderated F-statistic}
#'   \item{.p.value}{p-value generated from moderated F-statistic}
#'   \item{.df.total}{total degrees of freedom per gene}
#'   \item{.df.residual}{residual degrees of freedom per gene}
#'   \item{.s2.post}{posterior estimate of residual variance}
#'
#' @export
setMethod("augment", "MArrayLM", function(x, data, ...) {
    cols <- c("sigma", "F", "F.p.value", "AMean", "df.total", "df.residual",
              "s2.prior", "s2.post")
    newnames <- c("sigma", "statistic", "p.value", "AMean", "df.total",
                  "df.residual", "s2.post")
    lst <- unclass(x)[cols]
    names(lst) <- paste0(".", newnames)
    ret <- cbind(.gene=rownames(x$coefficients), as.data.frame(compact(lst)))

    if (!missing(data)) {
        ret <- cbind(as.data.frame(data), ret)
    }

    finish(ret)
})


#' @rdname limma_tidiers
#'
#' @return \code{glance} returns one row, containing
#'   \item{rank}{rank of design matrix}
#'   \item{df.prior}{empirical Bayesian prior degrees of freedom}
#'   \item{s2.prior}{empirical Bayesian prior residual standard deviation}
#'
#' @export
setMethod("glance", "MArrayLM", function(x, ...) {
    components <- unclass(x)[c("rank", "df.prior", "s2.prior")]
    as.data.frame(compact(components))
})
