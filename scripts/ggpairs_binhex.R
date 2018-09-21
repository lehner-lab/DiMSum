
#ggpairs_binhex
#
# GGpairs plot for all variables in input data table with correlation in upper triangle and 2d binned hexagons in lower triangle.
#
# input_dt: input data table (required)
# output_file: plot output path (required)
# width: plot width (default: 12)
# height: plot height (default: 12)
# bins: number of hexagon bins (default: 50)
#
# Returns: nothing
#
ggpairs_binhex <- function(
  input_dt, 
  output_file, 
  width = 12, 
  height = 12, 
  bins = 50){
  d <- ggpairs(input_dt,
    columns = 1:dim(input_dt)[2],
    upper=list(continuous = "cor"),
    lower="blank", xlab = "log10(variant count + 1)", ylab = "log10(variant count + 1)")
  for (x in 1:dim(input_dt)[2]){
    for (y in 1:dim(input_dt)[2]){
      if (y>x) {
        d <- putPlot(d, ggplot(input_dt, aes_string(x=colnames(input_dt)[x],y=colnames(input_dt)[y])) + stat_binhex(bins=bins), y,x)
      }
    }
  }
  suppressWarnings(suppressMessages(ggsave(output_file, d, width = width, height = height)))
}
