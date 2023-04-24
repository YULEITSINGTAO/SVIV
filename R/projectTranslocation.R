#' Project translocation to Cartesian coordinate system
#'
#' @param translocation_bed dataframe, paired translocation bed dataframe
#' Columns:
#' -Chr_1 and Pos_1: First break point chromosome ID and location.
#' -Chr_2 and Pos_2: Second breal point chromosome ID and location.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' translocation_bed <- data.frame(Chr_1 = paste0("chr", c(7, 8, 9, 9, 13, 15)),
#' Pos_1 = c(8663231, 70602300, 131457166, 33130549, 21746650, 42750778),
#' Chr_2 = paste0("chr", c(5, 1, 2, 6, 11, 2)),
#' Pos_2 = c(37709720, 91853200, 116376668, 43655549, 108585748, 214996194))
#'
#' projectTranslocation(translocation_bed)
#'
projectTranslocation <- function(translocation_bed){

    translocation_bed$Event <- paste0("event_", c(1:nrow(translocation_bed)))
    translocation_bed$Chr_1 <- factor(translocation_bed$Chr_1, levels = paste0("chr", c(1:22, "X", "Y")))
    translocation_bed$Chr_2 <- factor(translocation_bed$Chr_2, levels = paste0("chr", c(1:22, "X", "Y")))

    p <- ggplot2::ggplot(translocation_bed, aes(Pos_1, Pos_2)) + geom_point() + facet_grid(vars(Chr_1), vars(Chr_2)) +
        scale_x_continuous(labels  = scales::label_number(scale = 1e-6, suffix = "Mbp", accuracy = 1)) +
        scale_y_continuous(labels  = scales::label_number(scale = 1e-6, suffix = "Mbp", accuracy = 1)) + xlab("Break point 1") + ylab("Break point 2") # for the y axis label


    return(p)
}
