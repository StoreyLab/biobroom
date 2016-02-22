#' Tidying methods for a sva list
#'
#' These are methods for turning a sva list, from the sva package, into a tidy data frame.
#' \code{tidy} returns a data.frame of the estimated surrogate variables, \code{glance} returns a data.frame
#' of the posterior probabilities, and \code{glance} returns a
#' data.frame with only the number of surrogate variables.
#'
#' @param x sva list
#' @param data Original data
#' @param ... extra arguments (not used)
#'
#' @return All tidying methods return a \code{data.frame} without rownames.
#' The structure depends on the method chosen.
#'
#' @rdname sva_tidiers
#'
#' @return \code{augment} returns one row per gene. It always
#' contains the columns
#'   \item{pprob.gam}{Posterior probability each gene is affected by heterogeneity}
#'   \item{pprob.b}{Posterior probability each gene is affected by model}
#' @export
augment_sva <- function(x, data, ...) {
    df <- data.frame(pprob.b = x$pprob.b, pprob.gam = x$pprob.gam)
    if (!missing(data)) {
        df <- cbind(as.data.frame(data), df)
    }
    df <- df[, !duplicated(colnames(df))]
    finish(df)
}

#' @rdname sva_tidiers
#'
#' @param addVar add additional coefficients to the estimated surrogate variables
#'
#' @return \code{tidy} returns the estimate surrogate variables.
#' @export
tidy_sva <- function(x, addVar = NULL, ...) {
    df <- data.frame(sv = x$sv)
    colnames(df) <- paste0("sv", 1:dim(df)[2])
    if (!is.null(addVar)) {
        if (nrow(addVar) != nrow(df)) stop("number of rows must be the same as
                                            the number of columns in original data")
        df <- cbind(as.data.frame(addVar), df)
    }
    df <- df[, !duplicated(colnames(df))]
    finish(df)
}
#' @rdname sva_tidiers
#'
#' @return \code{glance} returns the estimate surrogate variables.
#' @export
glance_sva <- function(x, ...) {
    df <- data.frame(n.sv = x$n.sv)
    finish(df)
}
