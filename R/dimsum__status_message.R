
#' dimsum__status_message
#'
#' Write status message to standard output.
#'
#' @param ... objects to be output as character vector (required)
#'
#' @return Nothing
#' @export
dimsum__status_message <- function(
	...
	){
	cat(sprintf(...), sep='', file=stdout())
}
