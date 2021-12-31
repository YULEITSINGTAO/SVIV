##################### Exported utility functions

#' Extra Homo sapiens hg38 genome coordinates
#' @description
#' - [hg38coord] returns a tibble of all hg38 chromosome coordinates
#' - [hg38special] returns a tibble of all hg38 special region coordinates
#' @return a tibble coordinates in bed-like format
#' @export
#'
#' @examples
#' hg38coord()
hg38coord <- function(){
    hg38 <- readr::read_tsv(system.file(package = "VCFComparison", "inst","extdata", "hg38.bed"),
                            col_names = FALSE, show_col_types = FALSE, progress = FALSE)
    names(hg38) <- c("Chrom", "Start", "End", "Abs_start", "Abs_end")
    structure(hg38, class = c("VCFComparison_hg38coord", "spec_tbl_df", "tbl_df", "tbl", "data.frame"))
}


#' @exportS3Method
print.VCFComparison_hg38coord <- function(x) {
    cat(
        "Returning human hg38 genome coordinates\n",
        tblue("Chrom:"), "Chromosome numbers, 1-22, X, Y\n",
        tblue("Start:"), "Chromosome start position\n",
        tblue("End:"), "Chromosome end position\n",
        tblue("Abs_start:"), "Chromosome absolute cumulative start position\n",
        tblue("Abs_end:"), "Chromosome absolute cumulative end position\n"
    )
    print(tibble::as_tibble(x))
    x
}


#' @rdname hg38coord
#' @export
#' @examples
#' hg38specail()
hg38special <- function(){
    hg38_special <- readr::read_tsv(system.file(package = "VCFComparison", "inst","extdata", "hg38_specail_regions.bed"),
                                    col_names = FALSE, show_col_types = FALSE, progress = FALSE)
    names(hg38_special) <- c("Chrom", "Start", "End", "Region")
    structure(hg38_special, class = c("VCFComparison_hg38special", "spec_tbl_df", "tbl_df", "tbl", "data.frame"))
}

#' @exportS3Method
print.VCFComparison_hg38special <- function(x) {
    cat(
        "Returning human hg38 genome coordinates\n",
        tblue("Chrom:"), "Chromosome numbers, 1-22, X, Y\n",
        tblue("Start:"), "Chromosome start position\n",
        tblue("End:"), "Chromosome end position\n",
        tblue("Region:"), "Chromosome special region type, like centromere, telomere, Low Mappability, etc.\n"
    )
    print(tibble::as_tibble(x))
    x
}

#' Get or set options for VCFComparison package
#' @param option string, what option to get or set
#' @param value any R object, the value you want to set for the option.
#' @details If value is not provided, get the current value of this option, if
#' non-NULL value is provided, set the option with this value. If the value is not
#' provided and the option value is unset, return `FALSE`
#'
#' ### Current possible values and defaults
#' - verbose: FALSE
#' @return see details
#' @export
#'
#' @examples
#' VCFComparisonOption("verbose") # get current value
#' VCFComparisonOption("verbose", TRUE) # set it to a new value
#' VCFComparisonOption("verbose") # check the value again
VCFComparisonOption <- function(option = NULL, value = NULL) {
    if(is.null(option)) {
        cat(
            tyellow("Possible VCFComparison package options are:\n"),
            tblue("verbose:    bool, verbosity level, default `FALSE`\n"),
            tblue("color_cont: string, continuous color gradient for plots, use `setColorContinuous` to set\n"),
            tblue("color_dis : string, discrete color palette for plots, use `setColorDiscrete` to set\n")
        )
        return(invisible())
    }
    stopifnot(is.character(option) && length(option) == 1)
    option <- glue("VCFComparison.{option}")
    if(is.null(value)) return(getOption(option, FALSE))
    do.call(options, stats::setNames(list(value), option))
}



#' VCFComparison package plot colors
#' @description get VCFComparison package plots colors, or set custom colors for
#' both continuous and discrete colors.
#' @details [setColorContinuous] will set color palettes for continuous plot gradients,
#' [setColorDiscrete] set the colors to use for discrete color palettes.
#' @param plot_colors character vector, colors to use for gradients for discrete palettes.
#'
#' If the `plot_colors` is keyword "current", only display current color options but **DO NOT**
#' set new color options.
#'
#' If the value is `"default"`, for [setColorDiscrete] it uses *Set3*, for
#' [setColorContinuous] it uses *plasma* gradient.
#'
#' @param n integer, for [setColorDiscrete] only, how many unique colors to create based on
#' provided colors. See [grDevices::colorRampPalette] for details.
#' @return nothing will be returned, by the plot color options will be set, use
#' `VCFComparisonOption("color_cont")` and `VCFComparisonOption("color_dis")`
#' to check current values after using setting the colors.
#' @param muted bool, muted the print message on console?
#' @export
#'
#' @examples
#' setColorContinuous(c("red", "blue", "yellow"))
#' setColorDiscrete(c("red", "blue", "yellow"), 10)
setColorContinuous <- function(plot_colors = "current", muted = FALSE) {
    stopifnot(is.character(plot_colors))
    stopifnot(is.logical(muted) && length(muted) == 1)
    plasma <- c("#F0F921FF", "#CC4678FF", "#0D0887FF")
    plot_colors <- if(all(plot_colors %in% "default")) plasma
                   else if (all(plot_colors %in% "current"))  VCFComparisonOption("color_cont")
                   else                                plot_colors
    lapply(plot_colors, function(x) {
        if(stringr::str_starts(plot_colors, "#", negate = TRUE) && (!x %in% colors()))
            logErr("Custom colors must be hex value or one of color names in `colors()`")
    })
    if(!muted) {
        unlist(lapply(grDevices::colorRampPalette(plot_colors)(getOption("width", 80)), function(x) crayon::make_style(x)("|"))) %>%
            cat(tblue("Continuous colors:\n"),., "\n", sep = "")
    }
    VCFComparisonOption("color_cont", plot_colors)
    invisible(plot_colors)
}

#' @rdname setColorContinuous
setColorDiscrete <- function(plot_colors = "current", n = 8, muted = FALSE) {
    stopifnot(is.character(plot_colors))
    stopifnot(is.numeric(n) && length(n) == 1 && n > 0)
    stopifnot(is.logical(muted) && length(muted) == 1)
    n <- as.integer(n)
    # set3 <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
    #           "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
    #           "#CCEBC5", "#FFED6F")
    set2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
    plot_colors <- if(all(plot_colors %in% "default")) set2
                   else if (all(plot_colors %in% "current"))  VCFComparisonOption("color_dis")
                   else                                plot_colors
    lapply(plot_colors, function(x) {
        if(stringr::str_starts(plot_colors, "#", negate = TRUE) && (!x %in% colors()))
            logErr("Custom colors must be hex value or one of color names in `colors()`", parentFrame = 3)
    })
    plot_colors <- grDevices::colorRampPalette(plot_colors)(n)
    if(!muted) {
        unlist(lapply(plot_colors, function(x) crayon::make_style(x)(x))) %>%
            cat(tblue("Discrete colors:\n"),., "\n")
    }
    VCFComparisonOption("color_dis", plot_colors)
    invisible(plot_colors)
}
############### not exported internal utility functions

#' Find out the parent calling function name
#'
#' @param parent_level int, how many parent call stack levels to traces back,
#' number must bigger than 1, because 0 is `getParentFrame` itself.
#' @param char bool, convert the return to character? `FALSE` will be original `call`
#' object
#' @details This function must be used in another function
#' @examples
#' abc <- function(){
#'     getParentFrame(1)
#' }
#'
#' def <- function(){
#'     abc()
#' }
#' abc()
#' def()
#' getParentFrame()
getParentFrame <- function(parent_level = 1, char = TRUE) {
    if(parent_level < 1L) stop("parent_level must bigger than 1")
    call_stack <- sys.calls()
    if(length(call_stack) == 1) stop("This function must be called within another function, are you in global level?")
    parent_level <- if(length(call_stack) <= parent_level) length(call_stack) - 1 else parent_level
    call_name <- sys.calls()[[sys.nframe()- parent_level]][1]
    if(char) as.character(call_name) else call_name
}



