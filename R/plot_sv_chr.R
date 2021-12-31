


#' Plot SVs by chromosome
#'
#' @param sv
#' @param label_karyo
#' @param ref string, what genome to use, default is "hg38"
#'
#' @return
#' @export
#'
#' @examples
svChrPlot <- function(
    sv,
    ref = "hg38",
    label_karyo =TRUE
){
 if(notFalsy(spsUtil::checkNameSpace("ggplot2"))) logErr("Install xxx")

}
