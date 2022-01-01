


#' Plot SVs by chromosome
#'
#' @param sv
#' @param label_karyo
#' @param ref string, what genome to use, default is "hg38"
#'
#' @return
#' @export
#' @details
#' You must install `karyoploteR` package to use this function. By default, this
#' package is not installed.
#' @examples
svChrPlot <- function(
    sv,
    ref = "hg38",
    label_karyo =TRUE
){
    checkPkg("karyoploteR")

    chr <- "chr6"
    alpha <- 0.5
    type_colors <- VCFComparisonOption("color_dis")
    # calc samples
    Samples <- sv$Sample %>% unique()
    total_track <- length(Samples)
    # plot base
    pp <- getDefaultPlotParams(plot.type = 1)
    pp$leftmargin <- 0.2
    pp$data1height <- 30*length(Samples)
    kp <- plotKaryotype(genome="hg38", chromosomes = chr, plot.params = pp)
    kpAddCytobandLabels(kp, force.all=TRUE, srt=90)
    lapply(seq_along(Samples), function(Sample_n) {
        # set up
        at <- autotrack(current.track = Sample_n, total.tracks = total_track)
        r0 <- at$r0
        r1 <- at$r1
        Sample <- Samples[Sample_n]
        sv_sub <- filter(sv, Chrom == chr, Sample == !!Sample)
        Types <- sv_sub$Type %>% unique()
        kpAddLabels(kp, labels=Sample, r0=r0, r1 = r1)

        # start to plot
        lapply(Types, function(Type) {
            sv_sub <- dplyr::filter(sv_sub, Type == !!Type)
            if(nrow(sv_sub) < 1) return()
            type_color <- switch(
                Type,
                "DUP" = type_colors[1],
                "DEL" = type_colors[2],
                "INS" = type_colors[3],
                "BND" = type_colors[5],
                "INV" = type_colors[4]
            )
            if(Type != "BND") {
                # rgb-A correction
                type_color <- col2rgb(type_color) %>%
                    {setNames(as.list(.), c("red", "green", "blue"))} %>%
                    {.[['alpha']] <- 255*alpha; .[['maxColorValue']] <- 255; .} %>%
                    {do.call(rgb, .)}
                regions <- glue('{sv_sub$Chrom}:{sv_sub$Start}-{sv_sub$End}')
                kpPlotRegions(kp, data=regions, col=type_color, r0=r0, r1 = r1, avoid.overlapping = FALSE)
                return(NULL)
            }
            kpPoints(kp, chr = chr, x = sv_sub$Start, pch="x", cex = 1.25, col=type_color, y = (r0 + r1)/2)
            NULL
        })
    }) %>% invisible()

}
