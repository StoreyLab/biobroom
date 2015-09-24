unrowname <- function(x) {
    rownames(x) <- NULL
    x
}

compact <- function(x) Filter(Negate(is.null), x)


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
