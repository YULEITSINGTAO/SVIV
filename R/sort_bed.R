## This is a function to order the translocation
## We keep the smaller index chr -> larger indes chr

sort_trans_bed <- function(input_trans_bed){
  input_trans_bed <- input_trans_bed %>% filter(Chr1!=Chr2)
  chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY")
  for (i in 1:nrow(input_trans_bed)) {
    chr1_factor <-factor(input_trans_bed[i,"Chr1"], levels=chrOrder)
    chr2_factor <-factor(input_trans_bed[i,"Chr2"], levels=chrOrder)
    if (as.numeric(chr1_factor) > as.numeric(chr2_factor)){
      a_1 <- input_trans_bed[i,"Chr2"]
      a_2 <- input_trans_bed[i,"Pos2"]

      input_trans_bed[i,"Chr2"] <- input_trans_bed[i,"Chr1"]
      input_trans_bed[i,"Pos2"] <- input_trans_bed[i,"Pos1"]
      input_trans_bed[i,"Chr1"] <- a_1
      input_trans_bed[i,"Pos1"] <- a_2
    }
  }

  input_trans_bed <- input_trans_bed %>% distinct()
  return(input_trans_bed)
}
