### Venn_detail### 

Venn_detail_process <- function(ven){
  ven_result <- ven@result
  split_ven_result <- strsplit(tolower(ven_result$Subset), "_")
  ven_sample_count <- c()
  
  for (i in 1:length(split_ven_result)){
    num <- str_count(split_ven_result[[i]], "sample")
    ven_sample_count = c(ven_sample_count, sum(num))
    
  }
  ven_result$samples_count <- ven_sample_count
  
  return(ven_result)
  
}
