#' Union of intervals
#' @description Union operation on two bed dataframes.
#' @param bed_df_1 dataframe, bed format dataframe noting one type of the interval mutation.
#' Columns:
#' - Chr:  string, chromosome numbers example: chr1, chr2, etc
#' - Start:  numeric integers, SV starting position
#' - End:    numeric integers, SV ending p1osition
#'
#' @param bed_df_2 dataframe, bed format dataframe noting one type of the interval mutation.
#' Columns:
#' - Chr:  string, chromosome numbers example: chr1, chr2, etc
#' - Start:  numeric integers, SV starting position
#' - End:    numeric integers, SV ending p1osition
#'
#' @param verbose bool, to display verbose information?
#'
#' @return bed_df
#' @export
#'
#' @examples
#' bed_df_1 <- data.frame(Chr = c(paste0("chr", c(1,1,1))), Start=c(100, 200, 300), End=c(150, 250, 350))
#' bed_df_2 <- data.frame(Chr = c(paste0("chr", c(1,1,1))), Start=c(120, 220, 320), End=c(150, 250, 350))
#' unionBedDf(bed_df_1, bed_df_2, verbose = FALSE)
#'
unionBedDf <- function(bed_df_1, bed_df_2, verbose = FALSE){
    union_df_list <- list()
    chromosomes <- paste0("chr", c(1:22,"X","Y"))
    for (chr in chromosomes) {
        chr_bed_df <- rbind(bed_df_1 %>% dplyr::filter(Chr == chr),  bed_df_2 %>% dplyr::filter(Chr == chr))

        if (nrow(chr_bed_df)==0){
            union_df_intervals <- data.frame(Chr = NA, Start = NA, End = NA)
            if (verbose == TRUE){
            print(paste("There is no SVs in", chr))
            }
        }else{
            chr_bed_df_intervals <- intervals::Intervals(chr_bed_df %>% dplyr::select(Start, End))
            union_df_intervals <- as.data.frame(intervals::interval_union(chr_bed_df_intervals))
            union_df_intervals <- cbind(chr, union_df_intervals)
            colnames(union_df_intervals) <- c("Chr", "Start", "End")
        }

        union_df_list[[chr]] <- union_df_intervals

    }
    union_df <- na.omit(do.call(rbind, union_df_list))
    rownames(union_df) <- NULL
    return(union_df)
}

#' Union of bed intervals list
#' @description
#' @param bed_list_1 dataframe, bed format dataframe noting one type of the interval mutation
#' @param bed_list_2 dataframe, bed format dataframe noting one type of the interval mutation
#' @param verbose bool, print verbose information?
#' @return bed_df
#' @export
#' @examples
#'
#' bed_df_1 <- data.frame(Chr = c(paste0("chr", c(1,1,1))), Start=c(100, 200, 300), End=c(150, 250, 350))
#' bed_df_2 <- data.frame(Chr = c(paste0("chr", c(1,1,1))), Start=c(120, 220, 320), End=c(150, 250, 350))
#' bed_list_1 <- list(bed_df_1, bed_df_2)
#' bed_list_2 <- list(bed_df_2, bed_df_1)
#'
#' unionBedList(bed_list_1, bed_list_2, verbose = FALSE)
#'
unionBedList <- function(bed_list_1, bed_list_2, verbose = FALSE){

    bed_list <- list()

    for (i in names(bed_list_1)){
        bed_list[[i]] <- unionBedDf(bed_list_1[[i]], bed_list_2[[i]])

    }
    return(bed_list)
}

#' @rdname unionBedDf
#' @export
#' @param bed_list_1 dataframe, bed format dataframe noting one type of the interval mutation
#' @param bed_list_2 dataframe, bed format dataframe noting one type of the interval mutation
#'
'%U%' <- function(bed_list_1, bed_list_2){
    unionBedList(bed_list_1, bed_list_2)
}
