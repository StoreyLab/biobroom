#' Tidiers for return values from functions that aren't S3 objects
#'
#' This method handles the return values of functions that return lists
#' rather than S3 objects, such as \code{sva},
#' and therefore cannot be handled by S3 dispatch.
#'
#' @param x list object
#' @param ... extra arguments, passed to the tidying function
#'
#' @details Those tiders themselves are implemented as functions of the
#' form tidy_<function> that are not exported.
#'
#' @name list_tidiers
#'
#' @export
tidy.list <- function(x, ...) {
    if (all( c("sv", "pprob.gam", "pprob.b", "n.sv")
            %in% names(x))) {
        # returned from sva
        tidy_sva(x, ...)
    } else {
        stop("No tidying method recognized for this list")
    }
}


#' @export
augment.list <- function(x, ...) {
    if (all( c("sv", "pprob.gam", "pprob.b", "n.sv")
             %in% names(x))) {
        # returned from sva
        augment_sva(x, ...)
    } else {
        stop("No augmenting method recognized for this list")
    }
}

#' @rdname list_tidiers
#'
#' @export
glance.list <- function(x, ...) {
    if (all( c("sv", "pprob.gam", "pprob.b", "n.sv")
             %in% names(x))) {
        # returned from sva
        glance_sva(x, ...)
    } else {
        stop("No glance method recognized for this list")
    }
}
