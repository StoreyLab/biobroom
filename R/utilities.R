#' strip rownames from an object
#'
#' @param x a data frame
unrowname <- function(x) {
    rownames(x) <- NULL
    x
}


#' Remove NULL items in a vector or list
#'
#' @param x a vector or list
compact <- function(x) Filter(Negate(is.null), x)
