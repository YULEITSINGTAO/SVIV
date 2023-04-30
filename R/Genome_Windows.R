#' Genome Windows
#'
#' @param reference character, hg19 or hg38
#' @param tilewidth numeric, size of tile
#'
#' @return
# @export
#'
#' @examples
#' GenomeWindows("hg19", 500000)
#'
#' @noRd
GenomeWindows <- function(reference = c("hg19", "hg38"), tilewidth = 500000){

    chrSizes <- GenomeInfoDb::getChromInfoFromUCSC(reference) %>% dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y")))
    chrSizes <- chrSizes$size
    names(chrSizes) <- paste0("chr", c(1:22, "X", "Y"))
    bins <- GenomicRanges::tileGenome(chrSizes, tilewidth = tilewidth, cut.last.tile.in.chrom=T)
    bins <- as.data.frame(bins)[,c(1,2,3)]
    colnames(bins) <- c("Chr", "Start", "End")
    return(bins)
}
