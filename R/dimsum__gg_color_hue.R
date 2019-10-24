
#' dimsum__gg_color_hue
#'
#' ggplot2-like colour scale in HCL space.
#'
#' @param n an integer number of colours required (required)
#'
#' @return a vector of colours
#' @export
dimsum__gg_color_hue <- function(
  n
  ){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

