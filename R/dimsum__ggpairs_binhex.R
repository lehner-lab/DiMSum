
#' dimsum__ggpairs_binhex
#'
#' GGpairs plot for all variables in input data table with correlation in upper triangle and 2d binned hexagons in lower triangle.
#'
#' @param input_dt input data table (required)
#' @param output_file_prefix plot output path prefix (required)
#' @param width plot width (default: 12)
#' @param height plot height (default: 12)
#' @param bins number of hexagon bins (default: 50)
#' @param xlab the title of the x-axis
#' @param ylab the title of the y-axis
#' @param title the plot title
#' @param cut a factor variable to colour the hexagon outlines (optional)
#' @param size line width of hex bins when 'cut' specified (default: 0.5)
#' @param thresholds named list of thresholds where item is numeric (list) of thresholds and name is linetype (optional)
#'
#' @return Nothing
#' @export
dimsum__ggpairs_binhex <- function(
  input_dt, 
  output_file_prefix, 
  width = 12, 
  height = 12, 
  bins = 50,
  xlab = "x",
  ylab = "y",
  title = "",
  cut = NULL,
  size = 0.5,
  thresholds = NULL){
  #Check if something to plot
  if(dim(input_dt)[1] == 0){
    warning("dimsum__ggpairs_binhex.R: No data to plot (empty data.table 'input_dt').", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    return(NULL)
  }
  d <- GGally::ggpairs(input_dt,
    columns = 1:dim(input_dt)[2],
    upper=list(continuous = "cor"),
    lower="blank", xlab = xlab, ylab = ylab, title = title)
  if(is.null(cut)){
    plot_dt <- input_dt
  }else{
    plot_dt <- input_dt[, cut := cut]
  }
  for (x in 1:dim(input_dt)[2]){
    for (y in 1:dim(input_dt)[2]){
      if (y>x) {
        dpiece <- ggplot2::ggplot(plot_dt, ggplot2::aes_string(x = colnames(input_dt)[x], y = colnames(input_dt)[y]))
        if(!is.null(thresholds)){
          dpiece <- dpiece + 
            ggplot2::geom_hline(yintercept = unlist(thresholds), linetype = unlist(names(thresholds)), color = "darkgrey") + 
            ggplot2::geom_vline(xintercept = unlist(thresholds), linetype = unlist(names(thresholds)), color = "darkgrey")
        }
        if(is.null(cut)){
          dpiece <- dpiece + ggplot2::stat_binhex(bins = bins, size = size, color = "lightgrey")
        }else{
          dpiece <- dpiece + ggplot2::stat_binhex(ggplot2::aes(color = cut), bins = bins, size = size)
        }
        dpiece <- dpiece + 
          ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
          ggplot2::geom_abline(color = "darkgrey", lty = 2)
      d <- GGally::putPlot(d, dpiece, y, x)    
      }
    }
  }
  d <- d + ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  suppressWarnings(suppressMessages(dimsum__save_png(paste0(output_file_prefix, ".png"), d, width = width, height = height)))
  suppressWarnings(suppressMessages(ggplot2::ggsave(paste0(output_file_prefix, ".pdf"), d, width = width, height = height)))
}
