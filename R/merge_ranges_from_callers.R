#' Merge ranges from different callers
#'
#' @param Manta_list x
#' @param DELLY_list xx
#' @param GRIDSS_list x
#' @param LUMPY_list x
#'
#' @return x
# @export
#'
#' @examples
#' 1 + 1
#' @noRd
MergeRangesFromCallers <- function(Manta_list, DELLY_list, GRIDSS_list, LUMPY_list){
  caller_overlap_grl <- list()
  for (i in 1:length(Manta_list)) {
    ## Convert the dataframe to GRange
    if(nrow(Manta_list[[i]]) == 0){
      Manta_gr = GRanges()
    }else{
      Manta_gr <- makeGRangesFromDataFrame(Manta_list[[i]])
    }

    if(nrow(DELLY_list[[i]]) == 0){
      DELLY_gr = GRanges()
    }else{
      DELLY_gr <- makeGRangesFromDataFrame(DELLY_list[[i]])
    }

    if(nrow(GRIDSS_list[[i]]) == 0){
      GRIDSS_gr = GRanges()
    }else{
      GRIDSS_gr <- makeGRangesFromDataFrame(GRIDSS_list[[i]])
    }

    if(nrow(LUMPY_list[[i]]) == 0){
      LUMPY_gr = GRanges()
    }else{
      LUMPY_gr <- makeGRangesFromDataFrame(LUMPY_list[[i]])
    }

    ## Get the overlap
    Manta_DELLY_overleap <- GenomicRanges::intersect(Manta_gr, DELLY_gr)
    GRIDSS_Manta_overlap <- GenomicRanges::intersect(GRIDSS_gr, Manta_gr)
    GRIDSS_DELLY_overlap <- GenomicRanges::intersect(GRIDSS_gr, DELLY_gr)
    LUMPY_Manta_overlap <- GenomicRanges::intersect(LUMPY_gr, Manta_gr)
    LUMPYS_DELLY_overlap <- GenomicRanges::intersect(LUMPY_gr, DELLY_gr)
    ## Reduce the GRange
    caller_overlap_grl[[i]] <- GenomicRanges::reduce(c(Manta_DELLY_overleap, GRIDSS_Manta_overlap,
                                                       GRIDSS_DELLY_overlap, LUMPY_Manta_overlap,
                                                       LUMPYS_DELLY_overlap))
  }
  return(caller_overlap_grl)
}
