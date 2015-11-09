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
#' @param x \code{MArrayLM}, \code{MAList}, \code{Elist} object
#'
#' @return The output of tidying functions is always a data frame without
#' rownames.
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
#' @param intercept whether the \code{(Intercept)} term should be included
#' (default FALSE)
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
#' @S3method tidy MArrayLM
#' @export tidy.MArrayLM
tidy.MArrayLM <- function(x, intercept = FALSE, ...) {
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
    if (!intercept) {
        ret <- ret %>% dplyr::filter(term != "(Intercept)") %>%
            dplyr::mutate(term = droplevels(term))
    }

    finish(ret)
}


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
#' @S3method augment MArrayLM
#' @export augment.MArrayLM
augment.MArrayLM <- function(x, data, ...) {
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
}


#' @rdname limma_tidiers
#'
#' @return \code{glance} returns one row, containing
#'   \item{rank}{rank of design matrix}
#'   \item{df.prior}{empirical Bayesian prior degrees of freedom}
#'   \item{s2.prior}{empirical Bayesian prior residual standard deviation}
#'
#' @S3method glance MArrayLM
#' @export glance.MArrayLM
glance.MArrayLM <- function(x, ...) {
    components <- unclass(x)[c("rank", "df.prior", "s2.prior")]
    as.data.frame(compact(components))
}


#' Tidying method for a MA list
#' @rdname limma_tidiers
#'
#' @return \code{tidy} returns a data frame with one row per gene-sample
#' combination, with columns
#'   \item{gene}{gene name}
#'   \item{sample}{sample name (from column names)}
#'   \item{value}{expressions on log2 scale}
#' @S3method tidy MAList
#' @export tidy.MAList
tidy.MAList <- function(x, ...) {
    tidy_matrix(x$M)
}


#' Tidy an EList expression object
#'
#' @rdname limma_tidiers
#' @param addTargets Add sample level information. Default is FALSE.
#' @return \code{tidy} returns a data frame with one row per gene-sample
#' combination, with columns
#'   \item{gene}{gene name}
#'   \item{sample}{sample name (from column names)}
#'   \item{value}{expressions on log2 scale}
#'   \item{weight}{present if \code{weights} is set}
#'   \item{other columns}{if present and if \code{addTargets} is set}
#' @S3method tidy EList
#' @export tidy.EList
tidy.EList <- function(x, addTargets=FALSE, ...) {
  ret <- tidy_matrix(x$E)
  if (!is.null(x$weights)) {
    rownames(x$weights) <- rownames(x$E)
    ret$weight <- tidy_matrix(x$weights)$value
  }
  if (!is.null(x$sample.weights)) {
      sw <- setNames(x$sample.weights, colnames(x))
      ret$sample.weight <- sw[ret$sample]
  }
  if (addTargets) {
      targets <- x$targets
      rownames(targets) <- colnames(x)
      ret <- cbind(
          ret[, setdiff(names(ret), 'value')],
          targets[ret$sample,,drop=FALSE],
          value=ret$value) %>%
          unrowname
  }
  ret
}

tidy_matrix <- function(x, ...) {
    broom::fix_data_frame(x, newcol = "gene") %>%
        tidyr::gather(sample, value, -gene) %>%
        dplyr::mutate(sample=as.character(sample))
}
