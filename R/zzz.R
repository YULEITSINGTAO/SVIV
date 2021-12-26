
## when namespace is loaded do following:
.onLoad <- function(libname, pkgname) {
    #default options
    VCFComparisonOption("verbose", FALSE)
    setColorContinuous(muted = TRUE)
    setColorDiscrete(muted = TRUE)
}
