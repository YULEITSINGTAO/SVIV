# Project internal logging methods

###########
#' internal functions for logging at
#' info - message level
#' warn - warning
#' err  - error
#' @noRd
#' @details
#' 1. Logging methods are meant to be used inside anthoer function, direct use
#'    is not recommended.
#' 2. by default info level is muted unless `verbose = TRUE`
#' 3. default prefix is the parent function call name, or you can specify a name
#'    you want.
#' @examples
#' abc <- function(){
#'     logInfo("mylog infooo") # this is muted
#'     logInfo("mylog", TRUE) # this is not muted
#' }
#'
#' def <- function() {
#'     logWarn("my warning")
#'     logErr("error happened")
#'     1 + 1 # this following code is blocked by the error
#' }
#' abc()
#' def()
logInfo <- function(
    msg,
    verbose = SVIVOption("verbose"),
    prefix="default",
    parentFrame = 2
    ){
    call_name <- if (prefix == "default") getParentFrame(parentFrame) else prefix
    if(verbose) msg(msg, info_text = if(notFalsy(call_name)) glue('{call_name}-INFO') else "INFO")
}

logWarn <- function(msg, prefix="default", parentFrame = 2){
    call_name <- if (prefix == "default") getParentFrame(parentFrame) else prefix
    msg(msg, level = "warning", warning_text = if(notFalsy(call_name)) glue('{call_name}-WARNING') else "WARNING")
}


logErr <- function(msg, prefix="default", parentFrame = 2){
    call_name <- if (prefix == "default") getParentFrame(parentFrame) else prefix
    msg(msg, level = "error", error_text = if(notFalsy(call_name)) glue('{call_name}-ERROR') else "ERROR")
}

############

############ different text colors
tgreen <- function(x){crayon::green$bold(x)}
tblue <- function(x){crayon::blue$bold(x)}
tred <- function(x){crayon::red$bold(x)}
tyellow <- function(x){crayon::yellow$bold(x)}
torange <- function(x){crayon::make_style("darkorange1")$bold(x)}



