
#' sample_count_distributions_by_ntmut
#'
#' Plot sample sample count distributions split by number of nucleotide mutations (Nmut_nt) and synonymous/nonsynonymous mutations (Nmut_aa).
#'
#' @param input_dt input data table (required)
#' @param output_file plot output path (required)
#' @param width plot width (default: 12)
#' @param height plot height (default: 8)
#' @param title plot title (default: "")
#'
#' @return Nothing
#' @export
sample_count_distributions_by_ntmut <- function(
  input_dt,
  output_file,
  width = 12,
  height = 8,
  title = ""){
  plot_df <- reshape2::melt(as.data.frame(input_dt), id = c("Nmut_nt", "Nmut_aa"))
  plot_df[,'Synonymous_variant'] <- factor(plot_df[,'Nmut_aa']==0)
  plot_df[,'Number_substitutions'] <- factor(plot_df[,'Nmut_nt'])
  plot_df[,'value'] <- plot_df[,'value']+1
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(value, color = Number_substitutions, linetype = Synonymous_variant)) +
    ggplot2::geom_density() + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
    ggplot2::labs(x = "variant count (log scale)", y = "Density", title = title) +
    ggplot2::facet_grid(Number_substitutions ~ variable, scales = "free")
  ggplot2::ggsave(output_file, width = width, height = height)
}

