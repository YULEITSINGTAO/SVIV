#' AnnotateIntervals
#'
#' This is function is used to annotate the gene ranges
#'
#' @param intervals GRange file
#' @param txdb gene database
#'
#' @return Annotated GRange
# @export
#'
#' @examples
#'
#' annotateIntervals(GRange, hg19_txdb)
#' @noRd
annotateIntervals <- function(intervals, txdb){

  stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
  broads <- GenomicFeatures::genes(txdb)
  anno <- GenomicRanges::nearest(intervals, broads)
  mcols(intervals)$nearest_ENTREZ_id[!is.na(anno)] <- mcols(broads[anno[!is.na(anno)]])$gene_id
  return(intervals)
  }
