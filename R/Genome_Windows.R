#' Genome Windows
#'
#' @param reference
#' @param tilewidth
#'
#' @return
# @export
#'
#' @examples
#' Genome_Windows("hg19", 500000)
#'
Genome_Windows <- function(reference = "hg19", tilewidth = 500000){

    if (reference == "hg19"){
        chrSizes <- readRDS(system.file(package = "VCFComparison", "extdata", "chr_length.rds"))

    } else if (reference == "hg38"){
        chrSizes <- hg38coord()
        chrSizes <- chrSizes$End
        names(chrSizes) <- paste0("chr", c(1:22, "X", "Y"))
    }

    bins <- GenomicRanges::tileGenome(chrSizes, tilewidth = tilewidth, cut.last.tile.in.chrom=T)
    bins <- as.data.frame(bins)[,c(1,2,3)]
    colnames(bins) <- c("Chr", "Start", "End")
    return(bins)
}
