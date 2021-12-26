
#' validate input SV dataframe
#' @param special df
#' @param show_BND bool, check BND?
#' @return checked results in a list
#' @noRd
.validateSVdf <- function(sv, show_BND) {
    logInfo("Validating input SV")
    has_group <- has_bnd_end <- TRUE
    if(!inherits(sv, "data.frame")) logErr("Input SVs must be a data.frame")
    sv_names <- names(sv)
    if(!all( c("Chrom", "Start", "End", "Type", "Sample") %in% sv_names))
       logErr('Input SV must have "Chrom", "Start", "End", "Type", "Sample" columns')
    if(!inherits(sv$Chrom, c("factor", "character"))) logErr('Column "Chrom" must be `factor` or `character`')
    if(!inherits(sv$Sample, c("factor", "character"))) logErr('Column "Sample" must be `factor` or `character`')
    if(!inherits(sv$Start, "numeric")) logErr('Column "Start" must be `numeric`')
    if(!inherits(sv$End, "numeric")) logErr('Column "End" must be `numeric`')
    if(!inherits(sv$Type, c("factor", "character"))) logErr('Column "Type" must be `factor` or `character`')
    if(!all(sv$Type %in% c("BND", "DEL", "DUP", "INS", "INV"))) logErr('Column "Type" values must be one of "BND", "DEL", "DUP", "INS", "INV"')
    if(! "Group" %in% sv_names) {
        logWarn("column 'Group' is recommend but is missing")
        has_group <- FALSE
    } else if(!inherits(sv$Group, c("factor", "character"))) logErr('Column "group" must be `factor` or `character`')

    if(show_BND) {
        if(! "Bnd_end" %in% sv_names) {
            logWarn("column 'Bnd_end' is info is required to mark translocations but missing, skip BND arc will not be shown")
            has_bnd_end <- FALSE
        } else {
            if(!inherits(sv$Bnd_end, c("character"))) logErr('Column "Bnd_end" must be `character`')
            bnd_values <- filter(sv, Type == "BND") %>% pull(Bnd_end)
            if(any(is.na(bnd_values))) logErr("Row value of 'Bnd_end' couln of BND entries cannot be `NA`")
            if(!all(stringr::str_detect(bnd_values, "\\w+[0+9]{0,2}:[0-9]+"))) logErr(
                'Values  in "Bnd_end" column of BND entries must be in format of "chromsome:poisition_number", e.g. "chr3:10000"'
            )
        }
    } else {has_bnd_end <- FALSE}


    logInfo("Validating input SV success")
    list(
        group = has_group,
        bnd_end = has_bnd_end
    )
}

#' validate input reference gnome coordinates
#' @param special df
#' @return checked df
#' @noRd
.validateRefCoord <- function(ref) {
    logInfo("Validating input reference genome coordinates")
    if(!inherits(ref, "data.frame")) logErr("Input reference must be a data.frame, use `hg38coord()` to see an example", parentFrame = 3)
    ref_names <- names(ref)
    if(!all( c("Chrom", "Start", "End") %in% ref_names))
        logErr('Input reference must have "Chrom", "Start", "End" columns')
    if(!inherits(ref$Chrom, c("character"))) logErr('Column "Chrom" must be `character`')
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
    dplyr::select(ref, Chrom, Start, End, Abs_start, Abs_end)
}

#' validate input reference special region coordinates
#' @param special df
#' @return checked df
#' @noRd
.validateSpecialRegion <- function(special) {
    logInfo("Validating input gnome special region coordinates")
    if(!inherits(special, "data.frame")) logErr("Input special regions must be a data.frame, use `hg38special()` to see an example", parentFrame = 3)
    special_names <- names(special)
    if(!all( c("Chrom", "Start", "End", "Region") %in% special_names))
        logErr('Input special regions must have "Chrom", "Start", "End", "Region" columns')
    if(!inherits(special$Chrom, c("character"))) logErr('Column "Chrom" must be `character`')
    if(!inherits(special$Start, "numeric")) logErr('Column "Start" must be `numeric`')
    if(!inherits(special$End, "numeric")) logErr('Column "End" must be `numeric`')
    if(!inherits(special$Region, "character")) logErr('Column "Region" must be `character`')
    logInfo("Validating input special region success")
    special
}

#' @param sv validated sv table
#' @param ref validated referebce genome table
#' @noRd
.calcSVabsCoord <- function(sv, ref) {
    logInfo("Calculating SV absolute positions ...")
    tibble::add_column(sv, Abs_start = NA, Abs_end = NA, .after = "End") %>%
        dplyr::mutate(
            Abs_start = ref$Abs_start[match(Chrom, ref$Chrom)] + Start,
            Abs_end = ref$Abs_start[match(Chrom, ref$Chrom)] + End
        )
}


#' @noRd
#' @param sv sv table
#'
#' @return BND only table with end location processed, 3 new columns
#'
#' 1. Bnd_end_1, string, chrom number of end position
#' 2. Bnd_end_2, numeric, chrom position number of end position
#' 3. Bnd_end_abs, numeric, absolute position of end position
.findBND <- function(sv, ref){
    sv %>%
        mutate(Bnd_end = stringr::str_split(Bnd_end, ":")) %>%
        tidyr::unnest_wider(col = c("Bnd_end"), names_sep = "_") %>%
        mutate(
            Bnd_end_2 = as.numeric(Bnd_end_2),
            Bnd_end_abs = ref$Abs_start[match(Bnd_end_1, ref$Chrom)] + Bnd_end_2
        )
}

#' Plot the overview of SVs
#' @description Plot SVs across all chromosomes on a big segment plot
#' @param sv_table dataframe, a table with all SVs in a `bed`-like format.
#' Columns required:
#' - chr: string, chromosome numbers
#' - Start: numeric integers, SV starting position
#' - End: numeric integers, SV ending p1osition
#' - Type: string, type of the SV, like `"INS"`, `"DEL"`, `"INV"`, etc.
#' - Sample: string, sample IDs
#' - Group: string, sample group IDs
#' @param gnome_interval path, string, where did the intervaldefault value `"default"
#' @param gnome_gap
#' @param gnome_blacklist
#'
#' @return
#' @export
#'
#' @examples
svOverviewPlot <- function(
    sv,
    ref = hg38coord(TRUE),
    spe = hg38special(TRUE),
    show_BND = TRUE,
    title = "SV overview plot",
    title_hjust = 0.4,
    xlab = "Chromosome",
    ylab = "Sample",
    alpha = 1,
    sample_block = TRUE,
    group_block = TRUE,
    xend_expand = 1.005
){
    logInfo("Inputs validating...")
    sv_vd_res <- .validateSVdf(sv, show_BND)
    ref <- .validateRefCoord(ref)
    .validateSpecialRegion(spe)
    sv <- .calcSVabsCoord(sv, ref)
    sv_bnd <- filter(sv, Type == "BND")
    sv <- filter(sv, Type != "BND")

    logInfo("Start to create plot base...") # BND not plotted at this point
    p <- ggplot() +
        geom_segment(data = sv, aes(
            x = Sample, xend = Sample,
            y = Abs_start, yend = Abs_end,
            color = Type
        ),
        size = 5,
        alpha = alpha
        # arrow = arrow(angle = 90, length = unit(0.001, "cm"), type = "open")
        ) +
        geom_abline(intercept = ref$Abs_end, slope =0, size = 0.1, alpha = 0.7) +
        ylab(ylab) +
        xlab(xlab) +
        ggtitle(title) +
        coord_flip() +
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
            labels = ref$Chrom,
            expand = c(0, 0),
            limit = c(0, max(hg38$Abs_end) * xend_expand)
        )

    # plot BND here
    if(nrow(sv_bnd) > 0 && show_BND) {
        logInfo("Adding BND points to plot...")
        p <- p +
            geom_point(data = sv_bnd, aes(x = Sample, y = Abs_start), shape = 4, size = 0.3)

        # If end position of BND is provided, add arc
        # TODO not ideal to add arc to plot, use this code in other places
        if(sv_vd_res$bnd_end) {
            logInfo("Adding BND arcs to plot...")
            sv_bnd <- .findBND(sv_bnd, ref)
            p <- p +
                geom_point(data = sv_bnd, aes(x = Sample, y = Bnd_end_abs), shape = 4, size = 0.3) +
                geom_curve(data = sv_bnd, aes(x = Sample, y = Abs_start, xend = Sample, yend = Bnd_end_abs), curvature = -0.03)
        }
    }

    if(sv_vd_res$group && group_block) {
        logInfo("Adding group info to plot...")
        sv_sample_group <- sv %>% dplyr::group_by(Sample, Group) %>% dplyr::count() %>%
            dplyr::arrange(Group)
        sample_order <- sv_sample_group$Sample
        sv_sample_count <- sv_sample_group %>% dplyr::ungroup() %>% dplyr::count(Group) %>% dplyr::arrange(Group)
        vline_count <- cumsum(sv_sample_count$n)
        p <- p +
            scale_x_discrete(limits= sample_order) +
            geom_vline(xintercept = vline_count[-length(vline_count)] + 0.5, size = 0.8)
    }
    if(sample_block) {
        logInfo("Adding sample block info to plot...")
        p <- p + geom_vline(xintercept = seq_len(length(unique(sv$Sample))) + 0.5, size = 0.3, linetype = 2)
    }
    p
}
