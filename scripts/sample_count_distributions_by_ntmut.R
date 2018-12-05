
#sample_count_distributions_by_ntmut
#
# Plot sample sample count distributions split by number of nucleotide mutations (Nmut_nt) and synonymous/nonsynonymous mutations (Nmut_aa).
#
# input_dt: input data table (required)
# output_file: plot output path (required)
# width: plot width (default: 12)
# height: plot height (default: 8)
#
# Returns: nothing
#
sample_count_distributions_by_ntmut <- function(
  input_dt,
  output_file,
  width = 12,
  height = 8,
  title = ""){
  plot_df <- melt(as.data.frame(input_dt), id = c("Nmut_nt", "Nmut_aa"))
  plot_df[,'Synonymous_variant'] <- factor(plot_df[,'Nmut_aa']==0)
  plot_df[,'Number_substitutions'] <- factor(plot_df[,'Nmut_nt'])
  plot_df[,'value'] <- plot_df[,'value']+1
  d <- ggplot(plot_df, aes(value, color = Number_substitutions, linetype = Synonymous_variant)) +
    geom_density() + scale_x_log10() + theme_bw() +
    labs(x = "log10(variant count + 1)", y = "Density", title = title) +
    facet_grid(Number_substitutions ~ variable, scales = "free")
  ggsave(output_file, width = width, height = height)
}

