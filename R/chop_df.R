## Chop the range

#' Chop ranges
#'
#' @param bed_df
#'
#' @return chop_df_list
# @export
#'
#' @examples
#' function(bed_df)
#'
Chop_df <- function(bed_df){
    chromosomes <- paste0("chr", c(1:22,"X","Y"))
    chop_df_list <- list()
    for (chr in chromosomes) {
        chr_bed_df <- bed_df %>% dplyr::filter(Chr == chr)
        if (nrow(chr_bed_df)==0){
            choped_df_intervals <- data.frame(Chr = NA, Start = NA, End = NA)
            logWarn(paste("There is no SVs in ", chr))

        }else{
            chr_bed_df_intervals <- intervals::Intervals(chr_bed_df %>% dplyr::select(Start, End))
            chr_bed_df_union <- as.data.frame(intervals::interval_union(chr_bed_df_intervals))
            for (i in 1:nrow(chr_bed_df_union)) {

                chr_bed_df_in_interval_i <- chr_bed_df %>% dplyr::filter(Start >= chr_bed_df_union[i,1], End <= chr_bed_df_union[i,2])
                points_in_interval_i <- sort(unique(c(chr_bed_df_in_interval_i[,2], chr_bed_df_in_interval_i[,3])))
                choped_df_intervals <-c()

                for (j in c(1:(length(points_in_interval_i) - 1))){
                    choped_df_intervals <- rbind(choped_df_intervals, c(points_in_interval_i[j], points_in_interval_i[j+1]))
                }
                choped_df_intervals <- cbind(chr, choped_df_intervals)
                colnames(choped_df_intervals) <- c("Chr", "Start", "End")
            }
        }
        chop_df_list[[chr]] = choped_df_intervals
    }
    return(chop_df_list)
}
