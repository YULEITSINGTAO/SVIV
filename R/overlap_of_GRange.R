#' Title
#'
#' @param grl_1 GRange 1
#' @param grl_2 GRange 2
#'
#' @return The overlap of gtl_1 and grl_2
# @export
#'
#' @examples
#'
#' overlap_of_GRange(grl_1,grl_2)
#' @noRd
overlap_of_GRange <- function(grl_1,grl_2){
  sample_number <- length(grl_1)
  for (k in 1:sample_number){
    overlap_grl[[paste0("sample_", k)]] <- GenomicRanges::intersect(gtl_1[[k]], gtl_2[[k]])
  }
}
