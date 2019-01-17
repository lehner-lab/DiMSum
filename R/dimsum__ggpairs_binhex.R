
#' dimsum__ggpairs_binhex
#'
#' GGpairs plot for all variables in input data table with correlation in upper triangle and 2d binned hexagons in lower triangle.
#'
#' @param input_dt input data table (required)
#' @param output_file plot output path (required)
#' @param width plot width (default: 12)
#' @param height plot height (default: 12)
#' @param bins number of hexagon bins (default: 50)
#' @param xlab the title of the x-axis
#' @param ylab the title of the y-axis
#' @param title the plot title
#'
#' @return Nothing
#' @export
dimsum__ggpairs_binhex <- function(
  input_dt, 
  output_file, 
  width = 12, 
  height = 12, 
  bins = 50,
  xlab = "x",
  ylab = "y",
  title = ""){
  d <- GGally::ggpairs(input_dt,
    columns = 1:dim(input_dt)[2],
    upper=list(continuous = "cor"),
    lower="blank", xlab = xlab, ylab = ylab, title = title)
  for (x in 1:dim(input_dt)[2]){
    for (y in 1:dim(input_dt)[2]){
      if (y>x) {
        d <- GGally::putPlot(d, ggplot2::ggplot(input_dt, ggplot2::aes_string(x=colnames(input_dt)[x],y=colnames(input_dt)[y])) + ggplot2::stat_binhex(bins=bins), y,x)
      }
    }
  }
  d <- d + ggplot2::theme_bw()
  suppressWarnings(suppressMessages(ggplot2::ggsave(output_file, d, width = width, height = height)))
}
