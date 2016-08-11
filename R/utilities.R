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
# From broom package
fix_data_frame <- function(x, newnames = NULL, newcol = "term", check.names = FALSE)
{
    if (!is.null(newnames) && length(newnames) != ncol(x)) {
        stop("newnames must be NULL or have length equal to number of columns")
    }
    if (all(rownames(x) == seq_len(nrow(x)))) {
        ret <- data.frame(x, stringsAsFactors = FALSE)
        if (!is.null(newnames)) {
            colnames(ret) <- newnames
        }
    }
    else {
        ret <- data.frame(...new.col... = rownames(x), unrowname(x),
                          stringsAsFactors = FALSE, check.names = check.names)
        colnames(ret)[1] <- newcol
        if (!is.null(newnames)) {
            colnames(ret)[-1] <- newnames
        }
    }
    unrowname(ret)
}
