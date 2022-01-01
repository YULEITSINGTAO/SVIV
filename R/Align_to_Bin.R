#' Align_to_Bin
#'
#' @param bed_file
#' @param bin_reference
#'
#' @return
#' @export
#'
#' @examples
#' Align_to_Bin(bed_file, bin_reference)
Align_to_Bin <- function(bed_file, bin_reference){

    trans_DELLY_bin <- list()

    for (i in 1:length(DELLY_trans_list)) {
        trans_DELLY_bin[[i]] <- bins
        k <- 4
        for (j in 1:nrow(DELLY_trans_list[[i]])) {
            trans_DELLY_bin[[i]][, k] <- 0
            trans_DELLY_bin[[i]][which(bins$Chr == DELLY_trans_list[[i]][j, "Chr1"] & bins$Start <= DELLY_trans_list[[i]][j, "Pos1"] & bins$End >= DELLY_trans_list[[i]][j, "Pos1"]), k] <- 1
            trans_DELLY_bin[[i]][which(bins$Chr == DELLY_trans_list[[i]][j, "Chr2"] & bins$Start <= DELLY_trans_list[[i]][j, "Pos2"] & bins$End >= DELLY_trans_list[[i]][j, "Pos2"]), k] <- 1
            k <- k + 1
        }
        a <- trans_DELLY_bin[[i]][, 4:ncol(trans_DELLY_bin[[i]])]
        a <- t(unique(t(a)))
        trans_DELLY_bin[[i]] <- cbind(trans_DELLY_bin[[i]][,1:3],a)
    }

    trans_paried_DELLY <- character()
    k <- 1
    for (i in 1:length(trans_DELLY_bin)) {
        for (j in 4:ncol(trans_DELLY_bin[[i]])) {
            trans_paried_DELLY[k] <- paste(which(trans_DELLY_bin[[i]][,j]==1)[1], which(trans_DELLY_bin[[i]][,j]==1)[2] , sep = "-")
            k <- k+1
        }
    }

    sort(table(trans_paried_DELLY))


}
