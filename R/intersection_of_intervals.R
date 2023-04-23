#' Intersection of intervals
#' @description
#' @param bed_df_1 xxx
#' @param bed_df_2 xxx
#' @param verbose bool, print verbose information?
#' @return xxx
#' @export
#' @examples
#' 1+1
intersectionBedDf <- function(bed_df_1, bed_df_2, verbose = FALSE){
    intersection_df_list <- list()
    chromosomes <- paste0("chr", c(1:22,"X","Y"))
    for (chr in chromosomes) {
        chr_bed_df_1 <- bed_df_1 %>% dplyr::filter(Chr == chr)
        chr_bed_df_2 <- bed_df_2 %>% dplyr::filter(Chr == chr)

        if (nrow(chr_bed_df_1)*nrow(chr_bed_df_2)==0){
            intersection_df_intervals <- data.frame(Chr = NA, Start = NA, End = NA)
            if (verbose == TRUE){
                print(paste("There is no SVs in", chr))
            }
        }else{

            chr_bed_df_intervals_1 <- intervals::Intervals(chr_bed_df_1 %>% dplyr::select(Start, End))
            chr_bed_df_intervals_2 <- intervals::Intervals(chr_bed_df_2 %>% dplyr::select(Start, End))

            intersection_df_intervals <- as.data.frame(intervals::interval_intersection(chr_bed_df_intervals_1, chr_bed_df_intervals_2))

            if(nrow(intersection_df_intervals) == 0){
                if (verbose == TRUE){
                    print(paste("There is no overlapping SVs in", chr))
                }
            }else{
                intersection_df_intervals <- cbind(chr, intersection_df_intervals)
                colnames(intersection_df_intervals) <- c("Chr", "Start", "End")
            }
        }

        intersection_df_list[[chr]] <- intersection_df_intervals

    }
    intersection_df <- na.omit(do.call(rbind, intersection_df_list))
    rownames(intersection_df) <- NULL
    return(intersection_df)

}

intersection_bed_list <- function(bed_list_1, bed_list_2, verbose = FALSE){

    bed_list <- list()

    for (i in names(bed_list_1)){
        bed_list[[i]] <- intersectionBedDf(bed_list_1[[i]], bed_list_2[[i]])

    }
    return(bed_list)
}

#' @rdname intersectionBedDf
#' @export
'%I%' <- function(bed_list_1, bed_list_2){
    intersection_bed_list(bed_list_1, bed_list_2)
}
