
#' dimsum__hamming_distance
#'
#' Calculate hamming distance of two strings of the same length.
#'
#' @param x input string 1 (required)
#' @param y input string 2 (required)
#'
#' @return hamming distance
#' @export
dimsum__hamming_distance <- function(x, y){
  return(sum(strsplit(x, "")[[1]]!=strsplit(y, "")[[1]]))
}
