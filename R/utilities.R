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


#' Finish a tidied dataset and format for returning
#'
#' Each tidy and augment function calls this function before
#' returning. It ensures the object has no rownames, then turns the
#' return object into a tbl_df, tbl_dt, or whatever form is requested
#' using the global biobroom.return option (default tbl_df)
#'
#' @param x a data frame
#'
#' @importFrom data.table as.data.table
#'
#' @return a data-frame-like object, such as a tbl_df,
#' tbl_dt, or data.table, whose class is determined by the
#' biobroom.return option
finish <- function(x) {
    x <- unrowname(x)

    opt <- getOption("biobroom.return", default = "tbl_df")
    if (opt == "tbl_df") {
        dplyr::tbl_df(x)
    } else if (opt == "tbl_dt") {
        dplyr::tbl_dt(x)
    } else if (opt == "data.table") {
        data.table::as.data.table(tbl_dt(x))
    } else if (opt == "data.frame") {
        as.data.frame(x)
    } else {
        stop(paste("Invalid biobroom.return format", opt))
    }
}
