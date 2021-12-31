
## when namespace is loaded do following:
.onLoad <- function(libname, pkgname) {
    #default options
    if(!VCFComparisonOption("verbose")) VCFComparisonOption("verbose", FALSE)
    if(isFalsy(VCFComparisonOption("color_cont"))) setColorContinuous("default", muted = TRUE)
    if(isFalsy(VCFComparisonOption("color_dis"))) setColorDiscrete("default", muted = TRUE)
}
