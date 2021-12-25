SV_Venn <- function(gr_list){
  library(GenomicRanges)
  vector_source <- c(1:length(gr_list))
  overlap_gr_list <- GRangesList()
  for (k in 2:length(vector_source)){
    print(k)
    combination_matrix <- combn(vector_source, k, FUN = NULL, simplify = TRUE)
    for (i in 1:ncol(combination_matrix)){
      overlap_gr <- gr_list[[combination_matrix[1,i]]]
      for (j in 1:nrow(combination_matrix)){
        overlap_gr <- GenomicRanges::intersect(overlap_gr, gr_list[[combination_matrix[j,i]]])
        overlap_gr
      }
      if(length(overlap_gr)!=0){
        overlap_gr_list[[toString(combination_matrix[,i])]] <- overlap_gr
      }
    }
  }
  return(overlap_gr_list)
}
