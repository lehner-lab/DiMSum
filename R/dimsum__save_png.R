
#' dimsum__save_png
#'
#' Save plot object as PNG.
#'
#' @param output_file output file path (required)
#' @param plot_object plot object (required)
#' @param width plot width in inches (default:7)
#' @param height plot height in inches (default:7)
#' @param plot_dpi dots per inch (default:300)
#'
#' @return Nothing
#' @export
dimsum__save_png <- function(
  output_file, 
  plot_object, 
  width = 7, 
  height = 7,
  plot_dpi = 300){

  #Open device
  Cairo::CairoPNG(
    output_file, 
    width = width, 
    height = height, 
    units = "in", 
    dpi = plot_dpi,
    bg = "white")
  #Plot
  if(names(plot_object)[1]=="grobs"){
    plot(plot_object)
  }else{
    print(plot_object)
  }
  #Close device
  { sink("/dev/null")
    dev.off()
    sink()
  }
}
