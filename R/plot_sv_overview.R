
#' validate input SV dataframe
#' @param special df
#' @param show_BND bool, check BND?
#' @return checked results in a list
#' @noRd
.validateSVdf <- function(sv) {
    logInfo("Validating input SV")
    has_group <- has_bnd_end <- TRUE
    if(!inherits(sv, "data.frame")) logErr("Input SVs must be a data.frame")
    sv_names <- names(sv)
    if(!all( c("Chr", "Start", "End", "Type", "Sample") %in% sv_names))
       logErr('Input SV must have "Chr", "Start", "End", "Type", "Sample" columns')
    if(!inherits(sv$Chr, c("factor", "character"))) logErr('Column "Chr" must be `factor` or `character`')
    if(!inherits(sv$Sample, c("factor", "character"))) logErr('Column "Sample" must be `factor` or `character`')
    if(!inherits(sv$Start, "numeric")) logErr('Column "Start" must be `numeric`')
    if(!inherits(sv$End, "numeric")) logErr('Column "End" must be `numeric`')
    if(!inherits(sv$Type, c("factor", "character"))) logErr('Column "Type" must be `factor` or `character`')
    if(!all(sv$Type %in% c("BND", "DEL", "DUP", "INS", "INV"))) logErr('Column "Type" values must be one of "BND", "DEL", "DUP", "INS", "INV"')
    if(! "Group" %in% sv_names) {
        logWarn("column 'Group' is recommend but is missing")
        has_group <- FALSE
    } else if(!inherits(sv$Group, c("factor", "character"))) logErr('Column "group" must be `factor` or `character`')

    # if(show_BND) {
    #     if(! "Bnd_end" %in% sv_names) {
    #         logWarn("column 'Bnd_end' is info is required to mark translocations but missing, skip BND arc will not be shown")
    #         has_bnd_end <- FALSE
    #     } else {
    #         if(!inherits(sv$Bnd_end, c("character"))) logErr('Column "Bnd_end" must be `character`')
    #         bnd_values <- filter(sv, Type == "BND") %>% pull(Bnd_end)
    #         if(any(is.na(bnd_values))) logErr("Row value of 'Bnd_end' couln of BND entries cannot be `NA`")
    #         if(!all(stringr::str_detect(bnd_values, "\\w+[0+9]{0,2}:[0-9]+"))) logErr(
    #             'Values  in "Bnd_end" column of BND entries must be in format of "chromsome:poisition_number", e.g. "chr3:10000"'
    #         )
    #     }
    # } else {has_bnd_end <- FALSE}


    logInfo("Validating input SV success")
    list(
        group = has_group
        # bnd_end = has_bnd_end
    )
}

#' validate input reference gnome coordinates
#' @param special df
#' @return checked df
#' @noRd
.validateRefCoord <- function(ref) {
    logInfo("Validating input reference genome coordinates")
    if(!inherits(ref, "data.frame")) logErr("Input reference must be a data.frame, use `refcoord()` to see an example", parentFrame = 3)
    ref_names <- names(ref)
    if(!all( c("Chr", "Start", "End") %in% ref_names))
        logErr('Input reference must have "Chr", "Start", "End" columns')
    if(!inherits(ref$Chr, c("character"))) logErr('Column "Chr" must be `character`')
    if(!inherits(ref$Start, "numeric")) logErr('Column "Start" must be `numeric`')
    if(!inherits(ref$End, "numeric")) logErr('Column "End" must be `numeric`')

    if(! "Abs_end" %in% ref_names) {
        logWarn("column 'Abs_end' is recommend but is missing, calculating")
        ref <- ref %>% dplyr::mutate(Abs_end = cumsum(ref$End))
    } else if(!inherits(ref$Abs_end, c("numeric"))) logErr('Column "Abs_end" must be `numeric`')
    if(! "Abs_start" %in% ref_names) {
        logWarn("column 'Abs_start' is recommend but is missing, calculating")
        ref <- ref %>% dplyr::mutate(Abs_start = c(0, Abs_end[-nrow(ref)] + 1))
    } else if(!inherits(ref$Abs_start, c("numeric"))) logErr('Column "Abs_start" must be `numeric`')
    logInfo("Validating input reference genome success")
    dplyr::select(ref, Chr, Start, End, Abs_start, Abs_end)
}

#' validate input reference special region coordinates
#' @param special df
#' @return checked df
#' @noRd
.validateSpecialRegion <- function(special) {
    logInfo("Validating input gnome special region coordinates")
    if(!inherits(special, "data.frame")) logErr("Input special regions must be a data.frame, use `hg38special()` to see an example", parentFrame = 3)
    special_names <- names(special)
    if(!all( c("Chr", "Start", "End", "Region") %in% special_names))
        logErr('Input special regions must have "Chr", "Start", "End", "Region" columns')
    if(!inherits(special$Chr, c("character"))) logErr('Column "Chr" must be `character`')
    if(!inherits(special$Start, "numeric")) logErr('Column "Start" must be `numeric`')
    if(!inherits(special$End, "numeric")) logErr('Column "End" must be `numeric`')
    if(!inherits(special$Region, "character")) logErr('Column "Region" must be `character`')
    logInfo("Validating input special region success")
    TRUE
}

#' @param sv validated sv table
#' @param ref validated referebce genome table
#' @noRd
.calcSVabsCoord <- function(sv, ref) {
    logInfo("Calculating SV absolute positions ...")
    tibble::add_column(sv, Abs_start = NA, Abs_end = NA, .after = "End") %>%
        dplyr::mutate(
            Abs_start = ref$Abs_start[match(Chr, ref$Chr)] + Start,
            Abs_end = ref$Abs_start[match(Chr, ref$Chr)] + End
        )
}


#' @noRd
#' @param sv sv table
#'
#' @return BND only table with end location processed, 3 new columns
#'
#' 1. Bnd_end_1, string, Chr number of end position
#' 2. Bnd_end_2, numeric, Chr position number of end position
#' 3. Bnd_end_abs, numeric, absolute position of end position
.findBND <- function(sv, ref){
    sv %>%
        mutate(Bnd_end = stringr::str_split(Bnd_end, ":")) %>%
        tidyr::unnest_wider(col = c("Bnd_end"), names_sep = "_") %>%
        mutate(
            Bnd_end_2 = as.numeric(Bnd_end_2),
            Bnd_end_abs = ref$Abs_start[match(Bnd_end_1, ref$Chr)] + Bnd_end_2
        )
}

#' Plot the overview of SVs
#' @description Plot SVs across all chromosomes and across all samples on a big segment plot
#' @param sv dataframe, a table with all SVs in a `bed`-like format dataframe (tibble).
#' Columns:
#' - Chr:  string, chromosome numbers
#' - Start:  numeric integers, SV starting position
#' - End:    numeric integers, SV ending p1osition
#' - Type:   string, type of the SV, must be "BND", "DEL", "DUP", "INS", "INV"
#' - Sample: string, sample IDs
#' - Group:  string, optional, sample group IDs
#'
#' @param ref dataframe, reference gnome
#' Columns:
#' - Chr:     string, Chromosome numbers,like  chr1-chr22, chrX, chrY
#' - Start:     numeric, Chromosome start position
#' - End:       numeric, Chromosome end position
#' - Abs_start: numeric, optional, Chromosome absolute cumulative start position
#' - Abs_end:   numeric, optional, Chromosome absolute cumulative end position
#'
#' Use [hg38coord()] to get an example of hg38 reference
#'
#' @param show_BND bool, plot translocation break points on the plot? Default is
#' `FALSE`. Showing a lot of BND can make the plot messy. BND break points are
#' marked with `X` symbols.
#' @param title string, plot title
#' @param title_hjust numeric, plot horizontal adjustment, default is trying to
#' move the title to the center of the plot
#' @param xlab string, x axis label
#' @param ylab  string, y axis label
#' @param alpha numeric, transparency of symbol and segments, between 0 and 1
#' @param sample_block bool, use dashed lines to separate between each sample?
#' @param group_block bool, if `Group` column presents in the `sv` table, use
#' solid lines to separate between each sample group?
#' @param xend_expand numeric, a fraction to expand the x-axis. Sometimes the last
#' chromosome is cropped, try increase this value a little bit so the last chromosome
#' can be displayed.
#' @param plot_margin vector of numbers, length of 4, the plot margin, unit in cm. The default is only used
#' to adjust the space on the left when group label is added, and it can prevent
#' overcrowding. If there is no group information, ignore this.
#' @param group_hjust numeric, a number to adjust the space between sample labels and
#' group labels
#' @param group_vjust numeric, a number to adjust the vertical alignment on Y-axis
#' @param color_palette string, color names or color hex codes of what colors you want
#' to use to mark SV types. The default is to follow package discrete color options.
#'
#' @return returns a ggplot object
#' @export
#'
#' @examples
#'
#' # read in simulated SV data
#' sv <- readr::read_csv(system.file(package = "SVIV","extdata", "Sim_SV.csv"))
#' # create a table with 3 patients and each with 2 samples
#' # for each sample, we give them 20 random SV events
#' sample_info <- tibble::tibble(
#'     Sample = rep(paste0("sample", 1:6), each = 20),
#'     Group = rep(paste0("patient", 1:3), each = 2 * 20)
#' )
#' # bind the SV and sample information together
#' set.seed(99)
#' sv <- sv %>%
#'     dplyr::slice_sample(n = 120) %>%
#'     # the size of SVs are too small to see on overview plot,
#'     # randomly add 1MB - 20MB length to each event start and end
#'     # so we can see it clearly for the eaxmple.
#'     # Please notice that in real life,
#'     # the size of SV varies and you do not typically see them
#'     # all in the MB-size range.
#'     dplyr::rowwise() %>%
#'     dplyr::mutate(
#'         Start = Start - as.integer(runif(1, 1e6, 2e7)),
#'         End = End - as.integer(runif(1, 1e6, 2e7))
#'     ) %>%
#'     dplyr::bind_cols(sample_info)
#'
#'
#' svOverviewPlot(sv)
svOverviewPlot <- function(
    sv,
    ref = hg38coord(),
    show_BND = FALSE,
    title = "SV overview plot",
    title_hjust = 0.4,
    ylab = "Chromosome",
    xlab = "",
    alpha = 1,
    sample_block = TRUE,
    group_block = TRUE,
    xend_expand = 1.005,
    plot_margin = c(0,0,0,0.5),
    group_hjust = -5,
    group_vjust = -0.5,
    color_palette = SVIVOption("color_dis")
){
    logInfo("Inputs validating...")
    stopifnot(is.logical(show_BND) && length(show_BND) == 1)
    stopifnot(is.character(title) && length(title) == 1)
    stopifnot(is.numeric(title_hjust) && length(title_hjust) == 1)
    stopifnot(is.character(xlab) && length(xlab) == 1)
    stopifnot(is.character(ylab) && length(ylab) == 1)
    stopifnot(is.numeric(alpha) && length(alpha) == 1)
    stopifnot(is.logical(sample_block) && length(sample_block) == 1)
    stopifnot(is.logical(group_block) && length(group_block) == 1)
    stopifnot(is.numeric(xend_expand) && length(xend_expand) == 1)
    stopifnot(is.character(color_palette))

    sv_vd_res <- .validateSVdf(sv)
    ref <- .validateRefCoord(ref)
    # spe_vd_res <- if(inherits(try(.validateSpecialRegion(spe)), "logical")) TRUE else FALSE
    # if(!spe_vd_res) logWarn("Special region validation failed, will not be plotted")
    sv <- .calcSVabsCoord(sv, ref)
    if(length(color_palette) < length(unique(sv$Type)))
        logErr("`color_palette` color values must be larger than types of SVs")

    sv_bnd <- filter(sv, Type == "BND")
    show_BND <- all(nrow(sv_bnd) > 0 && show_BND)
    sv <- filter(sv, Type != "BND")
    logInfo("Find out how many colors needed")
    p_colors <- color_palette[seq(unique(sv$Type) %>% length())]
    p_colors_n <- length(p_colors)
    names(p_colors) <- unique(sv$Type)
    if(show_BND) {
        p_colors <- c(p_colors, "BND" = color_palette[p_colors_n + 1])
        spsUtil::inc(p_colors_n)
    }

    logInfo("Start to create plot base...") # BND not plotted at this point
    p <- ggplot() +
        geom_abline(intercept = ref$Abs_end, slope =0, size = 0.1, alpha = 0.7) +
        ylab(ylab) +
        xlab(xlab) +
        ggtitle(title) +
        coord_flip(clip = "off") +
        theme_minimal() +
        theme(
            legend.position="bottom",
            axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = title_hjust)
        ) +
        scale_y_continuous(
            breaks = (ref$Abs_start + ref$Abs_end)/2,
            labels = ref$Chr,
            expand = c(0, 0),
            limits = c(0, max(ref$Abs_end) * xend_expand)
        ) +
        scale_color_manual(
            name = if(show_BND) "" else "Type",
            values = if(show_BND) p_colors[-p_colors_n] else p_colors,
            guide = guide_legend(override.aes = list(alpha = 1))
        )

    logInfo("Adding SVs...")
    p <- p + geom_segment(data = sv, aes(
        x = Sample, xend = Sample,
        y = Abs_start, yend = Abs_end,
        color = Type
    ),
    size = 5,
    alpha = alpha
    # arrow = arrow(angle = 90, length = unit(0.001, "cm"), type = "open")
    )
    # plot BND here
    if(show_BND) {
        logInfo("Adding BND points to plot...")
        p <- p +
            geom_point(data = sv_bnd, aes(x = Sample, y = Abs_start,  shape = Type),
                       size = 0.3, color = p_colors[p_colors_n], alpha = alpha) +
            scale_shape_manual(
                name = "Type", labels="BND", values = c(BND = 4),
                guide = guide_legend(override.aes = list(
                    size =  2, stroke = 2, alpha = 1
                ))
            )

        # If end position of BND is provided, add arc
        # TODO not ideal to add arc to plot, use this code in other places
        # if(sv_vd_res$bnd_end) {
        #     logInfo("Adding BND arcs to plot...")
        #     sv_bnd <- .findBND(sv_bnd, ref)
        #     p <- p +
        #         geom_point(data = sv_bnd, aes(x = Sample, y = Bnd_end_abs), shape = 4, size = 0.3) +
        #         geom_curve(data = sv_bnd, aes(x = Sample, y = Abs_start, xend = Sample, yend = Bnd_end_abs), curvature = -0.03)
        # }
    }

    if(sv_vd_res$group && group_block) {
        logInfo("Adding group info to plot...")
        sv_sample_group <- sv %>% dplyr::group_by(Sample, Group) %>% dplyr::count() %>%
            dplyr::arrange(Group)
        sample_order <- sv_sample_group$Sample
        sv_sample_count <- sv_sample_group %>% dplyr::ungroup() %>% dplyr::count(Group) %>% dplyr::arrange(Group)
        vline_count <- cumsum(sv_sample_count$n)
        # calculate where to add group label, the middle of samples in this group
        group_y_label = (vline_count + c(0, vline_count[-length(vline_count)]))/2

        p <- p +
            scale_x_discrete(limits= sample_order) +
            geom_vline(xintercept = vline_count[-length(vline_count)] + 0.5, size = 0.8) +
            ggplot2::geom_text(
                data = data.frame(),
                aes(y = 0, x = group_y_label, label = sv_sample_count$Group),
                size = 5, check_overlap = TRUE, angle = 90,
                hjust = group_vjust, vjust = group_hjust
            ) +
            theme(
                plot.margin = ggplot2::margin(plot_margin, unit = "cm")
            )
    }
    if(sample_block) {
        logInfo("Adding sample block info to plot...")
        p <- p + geom_vline(xintercept = seq_len(length(unique(sv$Sample))) + 0.5, size = 0.3, linetype = 2)
    }
    p
}
