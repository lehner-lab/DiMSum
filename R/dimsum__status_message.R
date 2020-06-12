
#' dimsum__status_message
#'
#' Write status message to standard output.
#'
#' @param ... objects to be output as character vector (required)
#' @param newline print newline after status message (default:F)
#'
#' @return Nothing
#' @export
dimsum__status_message <- function(
	...,
	newline=F
	){
	cat(sprintf(...), sep='', file=stdout())
	if(newline){
		cat("\n")
	}
}
