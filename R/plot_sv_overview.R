
#' validate input SV dataframe
#' @param special df
#' @return checked df
#' @noRd
.validateSVdf <- function(sv) {
    logInfo("Validating input SV")
    if(!inherits(sv, "data.frame")) logErr("Input SVs must be a data.frame")
    sv_names <- names(sv)
    if(!all( c("Chrom", "Start", "End", "Type", "Sample") %in% sv_names))
       logErr('Input SV must have "Chrom", "Start", "End", "Type", "Sample" columns')
    if(!inherits(sv$Chrom, c("factor", "character"))) logErr('Column "Chrom" must be `factor` or `character`')
    if(!inherits(sv$Sample, c("factor", "character"))) logErr('Column "Sample" must be `factor` or `character`')
    if(!inherits(sv$Start, "numeric")) logErr('Column "Start" must be `numeric`')
    if(!inherits(sv$End, "numeric")) logErr('Column "End" must be `numeric`')
    if(!inherits(sv$Type, c("factor", "character"))) logErr('Column "Type" must be `factor` or `character`')
    if(! "Group" %in% sv_names) logWarn("column 'Group' is recommend but is missing")
    else if(!inherits(sv$Group, c("factor", "character"))) logErr('Column "group" must be `factor` or `character`')
    logInfo("Validating input SV success")
    sv
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
.calcSVabsCoord <- function(sv, ref) {
    logInfo("Calculating SV absolute positions ...")
    tibble::add_column(sv, Abs_start = NA, Abs_end = NA, .after = "End") %>%
        dplyr::mutate(
            Abs_start = ref$Abs_start[match(Chrom, ref$Chrom)] + Start,
            Abs_end = ref$Abs_start[match(Chrom, ref$Chrom)] + End
        )
}

#' Plot the overview of SVs
#' @description Plot SVs across all chromosomes on a big segment plot
#' @param sv_table dataframe, a table with all SVs in a `bed`-like format.
#' Columns required:
#' - chr: string, chromosome numbers
#' - Start: numeric integers, SV starting position
#' - End: numeric integers, SV ending position
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
    alpha = 1

){
    logInfo("Inputs validating...")
    .validateSVdf(sv)
    ref <- .validateRefCoord(ref)
    .validateSpecialRegion(spe)
    sv <- .calcSVabsCoord(sv, ref)

    logInfo("start to plot...")
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
        ylab("Chromosome") +
        xlab("Samples") +
        coord_flip() +
        theme_minimal() +
        theme(
            legend.position="bottom",
            axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        scale_y_continuous(breaks = (ref$Abs_start + ref$Abs_end)/2, labels = ref$Chrom, expand = c(0,0))

    if("Group" %in% names(sv)) {
        logInfo("Adding group info...")
        sv_sample_group <- sv %>% dplyr::group_by(Sample, Group) %>% dplyr::count() %>%
            dplyr::arrange(Group)
        sample_order <- sv_sample_group$Sample
        sv_sample_count <- sv_sample_group %>% ungroup() %>% count(Group) %>% dplyr::arrange(Group)
        vline_count <- cumsum(sv_sample_count$n)
        p <- p +
            scale_x_discrete(limits= sample_order) +
            geom_vline(xintercept = vline_count + 0.5, size = 0.8)+
            geom_vline(xintercept = seq_len(length(unique(sv$Sample))) + 0.5, size = 0.3, linetype = 2)
    }
    p
}
