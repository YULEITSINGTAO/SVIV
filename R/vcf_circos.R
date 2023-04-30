#' Circos plot from VCF file
#'
#' @param SV_fix
#'
#' @return circos plot
# @export
#'
#' @examples
#' vcf_circos(VCF_file)
#' @noRd
#'
vcf_circos <- function(SV_fix){
  ############## Translation table ################
  circos.clear()
  Translocation_table <- SV_fix%>%filter(grepl("BND",ID))
  chr_end <- as.data.frame(str_extract(string = Translocation_table[,5], pattern = "(?<=\\[).*(?=\\[)|(?<=\\]).*(?=\\])")
                           %>%str_split_fixed(":", 2))

  chr_start <- Translocation_table[,1:2]

  chr_start[,2] <- as.numeric(chr_start[,2])
  chr_end[,2] <- as.numeric(chr_end[,2])

  translocation_start <- cbind(chr_start,chr_start[,2] + 5)
  translocation_end <- cbind(chr_end,chr_end[,2] + 5)
  ##################################################################################


  ############## Deletion table ################
  Deletion_table <- SV_fix %>% filter(ALT =="<DEL>")
  bed_Delection_table <- cbind(Deletion_table[,1:2],str_extract(string = Deletion_table[,8], pattern = "(?<=END=)\\d+"))

  names(bed_Delection_table) <- c("Chr", "Start", "End")
  bed_Delection_table$Start <- as.numeric(bed_Delection_table$Start)
  bed_Delection_table$End <- as.numeric(bed_Delection_table$End)
  bed_Delection_table <- bed_Delection_table %>% add_column(Value = 1)

  #############Duplication table######################################################
  Duplication_table <- SV_fix %>% filter(ALT =="<DUP>"| ALT =="<DUP:TANDEM>")
  bed_Duplication_table <- cbind(Duplication_table[,1:2],str_extract(string = Duplication_table[,8], pattern = "(?<=END=)\\d+"))

  names(bed_Duplication_table) <- c("Chr", "Start", "End")
  bed_Duplication_table$Start <- as.numeric(bed_Duplication_table$Start)
  bed_Duplication_table$End <- as.numeric(bed_Duplication_table$End)
  bed_Duplication_table <- bed_Duplication_table %>% add_column(Value = 2)
  ###############################################################################


  #############Inversion table#################################################


  Inversion_table <- SV_fix %>% filter(ALT =="<INV>")
  bed_Inversion_table <- cbind(Inversion_table[,1:2],str_extract(string = Inversion_table[,8], pattern = "(?<=END=)\\d+"))

  names(bed_Inversion_table) <- c("Chr", "Start", "End")
  bed_Inversion_table$Start <- as.numeric(bed_Inversion_table$Start)
  bed_Inversion_table$End <- as.numeric(bed_Inversion_table$End)
  bed_Inversion_table <- bed_Inversion_table %>% add_column(Value = 3)
  #############################################################################

  bed_list <- list(bed_Delection_table, bed_Duplication_table, bed_Inversion_table)
  #############################################################################
  circos.par("start.degree" = 90)
  circos.par("track.height" = 0.05)
  circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
  circos.initializeWithIdeogram(species = "hg19")

  color <- c("green", "blue", "red")
  circos.genomicTrack(bed_list, stack = TRUE, bg.border = NA, bg.col = "#EFEFEF",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = i+0.5, ybottom = i-0.5,
                                           col = color[i], border = NA, ...)
                      })

  for (i in 1:dim(translocation_start)[1]) {
    circos.link(translocation_start[i,1],translocation_start[i,2],translocation_end[i,1],translocation_end[i,2])
  }


}
