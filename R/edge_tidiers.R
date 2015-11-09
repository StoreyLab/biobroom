#' Tidying methods for edge's deSet object
#'
#' @param x deSet object
#' @param data Original data can be added. Default is NULL.
#' @param addPheno whether columns should be included in the tidied output
#' for those in the ExpressionSet's phenoData
#' @param ... extra arguments (not used)
#'
#' @details \code{addPheno=TRUE} adds columns that are redundant (since they
#' add per-sample information to a per-sample-per-gene data frame), but that
#' are useful for some kinds of graphs and analyses.
#'
#' @rdname edge_tidiers
#'
#'
#' @return \code{tidy} returns a data frame with one row per gene-sample
#' combination, with columns
#'   \item{gene}{gene name}
#'   \item{sample}{sample name (from column names)}
#'   \item{value}{expressions on log2 scale}
#'
#' @S3method tidy deSet
#' @export tidy.deSet
tidy.deSet <- function(x, addPheno=FALSE, ...) {
  if (is.null(rownames(exprs(x))) | all(rownames(exprs(x)) == 1:nrow(exprs(x)))) {
    rownames(exprs(x)) <- paste0("g", 1:nrow(exprs(x)))
  }
  if (is.null(colnames(exprs(x))) | all(colnames(exprs(x)) == 1:ncol(exprs(x)))) {
    colnames(exprs(x)) <- paste0("sample", 1:ncol(exprs(x)))
  }
  expressions <- fix_data_frame(Biobase::exprs(x), newcol="gene")
  ret <- expressions %>%
      tidyr::gather(sample, value, -gene)%>%
      dplyr::mutate(sample=as.character(sample))


  if (addPheno) {
    pdat <- pData(x)
    rownames(pdat) <- colnames(exprs(x))
    ret <- unrowname(as.data.frame(cbind(gene=ret$gene,
                                         sample=ret$sample,
                                         pdat[ret$sample, , drop=FALSE],
                                         value=ret$value)))
  }
  finish(ret)
}



#' @rdname edge_tidiers
#'
#' @return \code{augment} returns a data.frame with
#'     \item{p.value}{the original p-values given to \code{qvalue}}
#'     \item{q.value}{the computed q-values}
#'     \item{lfdr}{the local false discovery rate}
#' @S3method augment deSet
#' @export augment.deSet
augment.deSet <- function(x, data, ...) {
    x <- x@qvalueObj
    if (is.null(x)) stop("qvalueObj slot is empty")
    augment(x)
}


#' @rdname edge_tidiers
#'
#' @return \code{glance} returns a data.frame with the model fits
#' @S3method glance deSet
#' @export glance.deSet
glance.deSet <- function(x, ...) {
    df <- data.frame(full.model = as.character(x@full.model)[2],
                     null.model = as.character(x@null.model)[2],
                     stringsAsFactors=FALSE)
    finish(df)
}

