#' Extract gene names from translocation events (merging translocation by genes)
#'
#' @param input_vcf_list vcf GRane list
#' @param txdb gene database
#'
#' @return translocation gene list
# @export
#'
#' @examples
#' translocation_gene_list(vcf_list, txdb)
#' @noRd
translocation_gene_list <- function(input_vcf_list, txdb){
  translocation_gene_list <- list()
  for (i in 1:length(input_vcf_list)){
    Translocation_table <- input_vcf_list[[i]]%>%filter(grepl("SVTYPE=BND",INFO))
    chr_end <- as.data.frame(str_extract(string = Translocation_table[,5], pattern = "(?<=\\[).*(?=\\[)|(?<=\\]).*(?=\\])")
                             %>%str_split_fixed(":", 2))
    chr_start <- Translocation_table[,1:2]
    chr_start[,2] <- as.numeric(chr_start[,2])
    chr_start[,3] <- chr_start[,2]
    names(chr_start) <- c("Chr","Start","End")
    chr_end[,2] <- as.numeric(chr_end[,2])
    chr_end[,3] <- chr_end[,2]
    names(chr_end) <- c("Chr","Start","End")
    chr_start <- chr_start %>% filter(grepl("chr[0-9]{1,}", Chr))
    chr_end <- chr_end %>% filter(grepl("chr[0-9]{1,}", Chr))
    Gr_start <- GenomicRanges::makeGRangesFromDataFrame(chr_start)
    Gr_end <- GenomicRanges::makeGRangesFromDataFrame(chr_end)
    Gr_start_with_gene <- annotateIntervals(Gr_start, txdb)
    Gr_end_with_gene <- annotateIntervals(Gr_end, txdb)
    translocation_gene_id <- c(mcols(Gr_start_with_gene)$nearest_gene_id, mcols(Gr_end_with_gene)$nearest_gene_id)
    translocation_gene_list[[paste0("sample-", i)]] <- unique(translocation_gene_id[!is.na(translocation_gene_id)])

  }
  return(translocation_gene_list)
}


