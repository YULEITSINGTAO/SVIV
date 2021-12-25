cat("Welcome to VCFComparision package development\n")
(function(){
    if(!interactive()) return(cat("Develop this package in an interactive session\n"))
    if(Sys.getenv("RSTUDIO", "") != "1") return(cat("We recommend to develop this package in Rstudio\n"))
    if(length(find.package("devtools", quiet = TRUE)) == 0) return(cat("package {devtools} is required but not installed\n"))
    if(length(find.package("roxygen2", quiet = TRUE)) == 0) return(cat("package {roxygen2} is required but not installed\n"))
    if(length(find.package("crayon", quiet = TRUE)) == 0) return(cat("package {crayon} is required but not installed\n"))
    if(length(find.package("pkgdown", quiet = TRUE)) == 0) cat("package {pkgdown} is recommended but not installed\n")
    devtools::load_all(".")
    cat(crayon::green$bold("All functions loaded\n"))
    cat(crayon::blue$bold("When functions are modified, use"),
        "Ctrl+Shift+L",
        crayon::blue$bold("to reload all."),
        "\n")
    cat(crayon::blue$bold("When roxygen help notes are written, use"),
        "Ctrl+Shift+D",
        crayon::blue$bold("to re-oxygenize."),
        "\n")
    cat(crayon::blue$bold("When development is done, use"),
        "Ctrl+Shift+E",
        crayon::blue$bold("to run R CMD check."),
        "\n")
})()
