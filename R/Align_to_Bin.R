#' AlignToBin
#' @description This function can count the range SVs in given ranges. The given ranges can be bins.
#' Bins can be equally divided, or the result form chopped dataframe.
#' @param bed_file
#' @param bin_reference
#'
#' @return Counts in each bin
# @export
#'
#' @examples
#' \dontrun{
#' Align_to_Bin(bed_file, bin_reference)
#'}
#'@noRd
#'
AlignToBin <- function(bed_file, bin_reference){

    B <- bin_reference

    B$count <- 0

    for (i in 1:nrow(bed_file)) {

        B[which((as.character(B[,1]) == bed_file[i,1]) & (B[,2] <= bed_file[i,3]) &
                  (B[,3] >= bed_file[i,2])), 4] <- B[which((as.character(B[,1]) == bed_file[i,1]) & (B[,2] <= bed_file[i,3]) &
                                                               (B[,3] >= bed_file[i,2])), 4] + 1
       }

    return(B)

}
