#'Tidying methods for GRanges and GRangesList objects.
#'
#' @param x GRanges or GRangesList object
#' @param ... Not used.
#'
#'@return All tidying methods return a \code{data.frame} without rownames. \code{tidy} returns one row for each range, which contains
#'\itemize{
#'\item{start} of the range
#'\item{end} of the range
#'\item{width} (or length) of the range
#'\item{names} of the range
#'\item{strand}
#'\item{seqname} Name of the sequence from which the range comes (usually the chromosome)
#'\item{metadata} Any included metadata, (ie, score, GC content)
#'}
#'For \code{GRangesList}, there will also be a column representing which group the ranges comes from.
#'\code{glance} returns a \code{data.frame} with the number of ranges, the number of sequences, and the number of groups (if applicable).
#'
#'@examples
#' if (require("GenomicRanges", "airway")) {
#' data(airway)
#'
#' # GRangesList object
#' air_gr <- rowRanges(airway)
#'
#' tidy(air_gr)
#' glance(air_gr)
#'
#' # GRanges object
#' air_gr <- rowRanges(airway)@@unlistData
#'
#' tidy(air_gr)
#' glance(air_gr)
#'
#'
#' }

#' @rdname GRanges_tidiers
#' @S3method tidy GRanges
#' @export tidy.GRanges
tidy.GRanges <- function(x, ...) {
  x.dt = as.data.frame(x@ranges)

  x.dt$strand = unlist(mapply(x@strand@values, x@strand@lengths, FUN = function(value, length) {
    as.character(rep(value, length))
  }), use.names=FALSE, recursive=FALSE)

  x.dt$seqname = unlist(mapply(x@seqnames@values, x@seqnames@lengths, FUN = function(value, length) {
    as.character(rep(value, length))
  }), use.names=FALSE, recursive=FALSE)

  x.dt = as.data.frame(cbind(x.dt, x@elementMetadata))

  finish(x.dt)
}


#' @rdname GRanges_tidiers
#' @S3method tidy GRangesList
#' @export tidy.GRangesList
tidy.GRangesList <- function(x, ...) {
  x.dt = tidy(x@unlistData)

  part.dt = as.data.frame(x@partitioning)

  x.dt$item = unlist(mapply(as.factor(part.dt$names), as.numeric(part.dt$width), SIMPLIFY = FALSE, FUN = function(value, length) {
    as.character(rep(value, length))
  }), use.names=FALSE, recursive=FALSE)
  finish(x.dt)
}


#' @rdname GRanges_tidiers
#' @S3method glance GRanges
#' @export glance.GRanges
glance.GRanges <- function(x, ...) {
  x.dt = tidy(x)
  ret = data.frame(ranges = length(gr),
                   sequences = length(unique(x.dt$seqname)))
  finish(ret)
}


#' @rdname GRanges_tidiers
#' @S3method glance GRangesList
#' @export glance.GRangesList
glance.GRangesList <- function(x, ...) {
  x.dt = tidy(x)
  ret = data.frame(ranges = length(gr),
                   sequences = length(unique(x.dt$seqname)),
                   lists = length(unique(x.dt$item)))
  finish(ret)
}
