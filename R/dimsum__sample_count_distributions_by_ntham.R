
#' dimsum__sample_count_distributions_by_ntham
#'
#' Plot sample sample count distributions split by nucleotide hamming distance (Nham_nt) and synonymous/nonsynonymous mutations (Nham_aa).
#'
#' @param input_dt input data table (required)
#' @param sequence_type coding potential of sequence: noncoding/coding (required)
#' @param output_file plot output path (required)
#' @param width plot width (default: 12)
#' @param height plot height (default: 8)
#' @param title plot title (default: "")
#'
#' @return Nothing
#' @export
dimsum__sample_count_distributions_by_ntham <- function(
  input_dt,
  sequence_type,
  output_file,
  width = 12,
  height = 8,
  title = ""){
  #Check if something to plot
  if(dim(input_dt)[1]==0){
    warning("dimsum__sample_count_distributions_by_ntham.R: No data to plot (empty data.table 'input_dt').", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    return(NULL)
  }
  plot_df <- reshape2::melt(as.data.frame(input_dt), id = c("Nham_nt", "Nham_aa"))
  plot_df[,'Synonymous_variant'] <- factor(plot_df[,'Nham_aa']==0)
  plot_df[,'Hamming_distance'] <- factor(plot_df[,'Nham_nt'])
  plot_df[,'value'] <- plot_df[,'value']+1
  if(sequence_type=="coding"){
    d <- ggplot2::ggplot(plot_df, ggplot2::aes(value, color = Hamming_distance, linetype = Synonymous_variant))
  }else{
    d <- ggplot2::ggplot(plot_df, ggplot2::aes(value, color = Hamming_distance))
  }
  d <- d + ggplot2::geom_density() + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
    ggplot2::labs(x = "variant count (log scale)", y = "Density", title = title) +
    ggplot2::facet_grid(Hamming_distance ~ variable, scales = "free")
  ggplot2::ggsave(output_file, width = width, height = height)
}

