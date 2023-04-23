
## when namespace is loaded do following:
.onLoad <- function(libname, pkgname) {
    #default options
    if(!SVIVOption("verbose")) SVIVOption("verbose", FALSE)
    if(isFalsy(SVIVOption("color_cont"))) setColorContinuous("default", muted = TRUE)
    if(isFalsy(SVIVOption("color_dis"))) setColorDiscrete("default", muted = TRUE)
}
