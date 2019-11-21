
#' dimsum__sample_count_distributions
#'
#' Plot sample nucleotide count distributions split by hamming distance (Nham_nt)
#'
#' @param input_dt input data table with count, WT and remaning column used for facetting (required)
#' @param output_file plot output path (required)
#' @param width plot width (default: 12)
#' @param height plot height (default: 8)
#' @param title plot title (default: "")
#' @param error_rate sequencing error rate (default: 1e-3)
#' @param seq_length sequence length 
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__sample_count_distributions <- function(
  input_dt,
  output_file,
  width = 12,
  height = 8,
  title = "",
  error_rate = 1e-3,
  seq_length){
  #Check if something to plot
  if(dim(input_dt)[1]==0){
    warning("dimsum__sample_count_distributions.R: No data to plot (empty data.table 'input_dt').", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    return(NULL)
  }

  #Column used for facetting
  facet_col <- names(input_dt)[!grepl("_count|WT", names(input_dt))]

  #Rename count columns
  names(input_dt)[grep("_count", names(input_dt))] <- dimsum__plot_samplename(names(input_dt)[grep("_count", names(input_dt))])

  #Real variants
  plot_dt <- reshape2::melt(input_dt, id = c(facet_col, "WT"))
  plot_dt[, Hamming_distance := factor(.SD[[1]]),,.SDcols = facet_col]
  plot_dt[,value := value+1]

  #WT frequency
  plot_dt_wt <- plot_dt[WT==T,.SD,,.SDcols = c("variable", "value")]

  #Fake single expected frequency
  plot_dt_fsingles <- plot_dt[WT==T]
  plot_dt_fsingles[, value := (round((value - 1)*error_rate/4) + 1)]
  plot_dt_fsingles[, Hamming_distance := factor(1)]

  #Fake double expected frequency distribution
  plot_dt_fdoubles <- plot_dt[Hamming_distance == 1]
  for(i in unique(plot_dt_fdoubles[,variable])){
    plot_dt_fdoubles[variable==i,value := (round(((value - 1)*(1 + seq_length*error_rate) - plot_dt_fsingles[variable==i,value - 1])*error_rate/4) + 1)]
  }
  plot_dt_fdoubles <- plot_dt_fdoubles[,.(value = median(value)),variable]
  plot_dt_fdoubles[, Hamming_distance := factor(2)]

  #Remove real variants with hamming distance ==0
  plot_dt <- plot_dt[get(facet_col)>0,]

  #Plot
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(value, color = Hamming_distance)) +
    ggplot2::geom_density() + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
    ggplot2::labs(x = "variant count + 1 (log scale)", y = "Density", title = title) +
    ggplot2::geom_vline(data = plot_dt_wt, ggplot2::aes(xintercept = value), linetype = 2) +
    ggplot2::geom_vline(data = plot_dt_fsingles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
    ggplot2::geom_vline(data = plot_dt_fdoubles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
    ggplot2::facet_grid(Hamming_distance ~ variable, scales = "free")
  ggplot2::ggsave(output_file, width = width, height = height)
}

