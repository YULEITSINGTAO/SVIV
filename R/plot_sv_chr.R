


#' Plot SVs by chromosome
#' @description Plot SVs at chromosome level.
#' @param sv dataframe, a table with all SVs in a `bed`-like format.
#' Columns:
#'
#' - Chr:  string, chromosome numbers
#' - Start:  numeric integers, SV starting position
#' - End:    numeric integers, SV ending p1osition
#' - Type:   string, type of the SV, must be "BND", "DEL", "DUP", "INS", "INV"
#' - Sample: string, sample IDs
#' - Group:  string, optional, sample group IDs
#'
#' @param label_karyo bool, label regions on the chromsome?
#' @param ref string, what genome to use, default is "hg38".
#' @param chr string, which chromosome to plot.
#' @param filter_sample character vector or `NULL`, what samples to plot instead of all samples.
#' @param filter_group character vector or `NULL`, what sample groups to plot instead of all groups.
#' This require you to have a `Group` column in `sv` table
#' @param mark_bnd bool, mark translocation as "X" on the plot? Usually there will
#' be a lot of BND events which may mess up the plot. Default is not to plot these
#' BNDs.
#' @param show_legend bool, show legend?
#' @param legend_x_mar numeric, what percentage value should the legend move vertically?
#' Default is to move downward below the karyotype. If you have many samples,
#' the legend may plotted with overlap of karyotype. Then you may want to change
#' this number.
#' @param alpha numeric, what is SV event block transparency value, between 0-1
#' @param plot_colors character vector, what plot colors to use. Default follows
#' package [setColorDiscrete] setting, see `?VCFComparisonOption`. You need to
#' provide at least 5 colors.
#' @param title string, plot title
#'
#' @return a `plotKaryotype` plot object, see `?plotKaryotype` for details
#' @export
#' @details
#' You must install `karyoploteR` package to use this function. By default, this
#' package is not installed.
#' @examples
#' # plot only when karyoploteR is installed
#' if(!length(spsUtil::checkNameSpace("karyoploteR", quietly = TRUE))) {
#'     # read in simulated SV data
#'     sv <- readr::read_csv(system.file(package = "VCFComparison","extdata", "Sim_SV.csv"))
#'     # create a table with 3 patients and each with 2 samples
#'     # for each sample, we give them 20 random SV events
#'     sample_info <- tibble::tibble(
#'         Sample = rep(paste0("Sample", 1:6), each = 5),
#'         Group = rep(paste0("Patient", 1:3), each = 2 * 5)
#'     )
#'     # bind the SV and sample information together
#'     set.seed(99)
#'     sv <- sv %>%
#'         dplyr::filter(Chr == "chr14") %>%
#'         dplyr::slice_sample(n = 30) %>%
#'         # the size of SVs are too small to see on overview plot,
#'         # randomly add 100KB - 1MB length to each event start and end
#'         # so we can see it clearly for the eaxmple.
#'         # Please notice that in real life,
#'         # the size of SV varies and you do not typically see them
#'         # all in the MB-size range.
#'         dplyr::rowwise() %>%
#'         dplyr::mutate(
#'             Start = Start - as.integer(runif(1, 1e5, 1e6)),
#'             End = End + as.integer(runif(1, 1e5, 1e6))
#'         ) %>%
#'         dplyr::bind_cols(sample_info)
#'
#'     svChrPlot(sv, chr = "chr14", title = "SVs on Chromosome 14")
#' }
svChrPlot <- function(
    sv,
    ref = "hg38",
    chr = "chr1",
    filter_sample = NULL,
    filter_group = NULL,
    label_karyo = TRUE,
    mark_bnd = FALSE,
    show_legend = TRUE,
    legend_x_mar = -0.1,
    alpha = 0.5,
    plot_colors = VCFComparisonOption("color_dis"),
    title = "SV plot by chromosome"

){
    # validations
    checkPkg("karyoploteR")
    # if(!isNamespaceLoaded("karyoploteR")) attachNamespace("karyoploteR")
    sv_has_group <- .validateSVdf(sv)$group
    logInfo("Filter by sample")
    if(!is.null(filter_sample)) {
        if(!is.character(filter_sample)) logErr("filter_sample must be a character vector.")
        if(!all(filter_sample %in% unique(sv$Sample))) {
            on.exit(cat(c(unique(sv$Sample), "\n"), sep = ", "), TRUE)
            logErr("Invalid filter_sample value, possible values are\n")
        }
        sv <- dplyr::filter(sv, Sample %in% filter_sample)
    }
    logInfo("Filter by group")
    if(!is.null(filter_group)) {
        if(!sv_has_group) logWarn("filter_group is provided but your input SV table doesn't have Group column, skip")
        else {
            if(!is.character(filter_group)) logErr("filter_group must be a character vector.")
            if(!all(filter_group %in% unique(sv$Group))) {
                on.exit(cat(c(unique(sv$Group), "\n"), sep = ", "), TRUE)
                logErr("Invalid filter_group value, possible values are\n")
            }
            sv <- dplyr::filter(sv, Group %in% filter_group)
        }
    }
    logInfo("Filter by chromosome")
    if(!is.character("chr") || length(chr) > 1) logErr("`chr` must be a length 1 character")
    if(!chr %in% paste0("chr", c(1:22, "X", "Y")))
        logWarn("Chromosome names are usually chr1-22, chrX or chrY. Ignore if you are not using human genome.")
    sv <- dplyr::filter(sv, chr == !!chr)
    logInfo("See if any sample left")
    if(nrow(sv) < 1) logErr("No sample left after filtering. Check your group and sample filter combination")
    logInfo("check plot colors")
    if(length(plot_colors) < 5) logErr("Your custom `plot_colors` must have at least 5 values")
    names(plot_colors)[seq(5)] <- c("DUP", "DEL", "INS", "INV", "BND")
    logInfo("other checks")
    stopifnot(is.logical(label_karyo) && length(label_karyo) == 1)
    stopifnot(is.logical(mark_bnd) && length(mark_bnd) == 1)
    stopifnot(is.logical(show_legend) && length(show_legend) == 1)
    stopifnot(is.numeric(legend_x_mar) && length(legend_x_mar) == 1)
    stopifnot(is.numeric(alpha) && length(alpha) == 1)
    stopifnot(is.character(title) && length(title) == 1)

    # calc samples
    Samples <- sv$Sample %>% unique()
    total_track <- length(Samples)
    # plot base
    pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
    pp$leftmargin <- 0.2
    pp$data1height <- 30*length(Samples)
    logInfo("Plot base")
    kp <- karyoploteR::plotKaryotype(genome=ref, chromosomes = chr, plot.params = pp, main = title)
    if(label_karyo) {logInfo("Add markers"); karyoploteR::kpAddCytobandLabels(kp, force.all=TRUE, srt=90)}
    if(show_legend) {
        logInfo("Add legends")
        total_events <- unique(sv$Type)
        graphics::legend(x = "bottom", fill = plot_colors[total_events],
               legend = names(plot_colors[total_events]),
               horiz=TRUE, inset=c(0, legend_x_mar),
               xpd=TRUE, box.lwd = 0, bg = "transparent")
    }
    # loop through each sample
    lapply(seq_along(Samples), function(Sample_n) {
        # set up
        at <- karyoploteR::autotrack(current.track = Sample_n, total.tracks = total_track)
        r0 <- at$r0
        r1 <- at$r1
        Sample <- Samples[Sample_n]
        Types <- sv$Type %>% unique()
        karyoploteR::kpAddLabels(kp, labels=Sample, r0=r0, r1 = r1)
        # start to plot
        lapply(Types, function(Type) {
            if(Type == "BND" && !mark_bnd) return(logInfo(c("BND skipped for ", Sample)))
            # get sample specific svs by type
            sv_sub <- dplyr::filter(sv, Type == !!Type & Sample == !!Sample)
            if(nrow(sv_sub) < 1) return(logInfo(c("No events of ", Type, " for ", Sample)))
            type_color <- plot_colors[Type]
            # for non-BND
            if(Type != "BND") {
                # rgb-A correction
                type_color <- ggplot2::alpha(type_color, alpha = alpha)
                regions <- glue('{sv_sub$Chr}:{sv_sub$Start}-{sv_sub$End}')
                karyoploteR::kpPlotRegions(kp, data=regions, col=type_color, r0=r0, r1 = r1, avoid.overlapping = FALSE)
                return(NULL)
            }
            # for BND
            karyoploteR::kpPoints(kp, chr = chr, x = sv_sub$Start, pch="x", cex = 1.25, col=type_color, y = (r0 + r1)/2)
            NULL
        })
    }) %>% invisible()
    kp
}
