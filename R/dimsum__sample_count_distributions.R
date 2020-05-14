
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
  fcols <- names(input_dt)[!grepl("_count|WT", names(input_dt))]

  #Add 1 + jitter
  set.seed(1)
  count_cols <- names(input_dt)[grepl("_count", names(input_dt))]
  for(i in count_cols){
    input_dt[, paste0(i) := as.numeric(.SD[[1]]+1),,.SDcols = i]
    input_dt[, paste0(i) := 10^(log10(.SD) + log10((.SD+1)/.SD) * runif(1)),1:nrow(input_dt) %% 1000,.SDcols = i]
    # input_dt[, paste0(i) := .SD + runif(1),1:nrow(input_dt) %% 1000,.SDcols = i]
  }

  #Rename count columns
  names(input_dt)[grep("_count", names(input_dt))] <- dimsum__plot_samplename(names(input_dt)[grep("_count", names(input_dt))])

  #Real variants
  plot_dt <- reshape2::melt(input_dt, id = c(fcols, "WT"))
  plot_dt[, Hamming_distance := factor(.SD[[1]]),,.SDcols = fcols]
  # plot_dt[,value := value+1]

  #WT frequency
  plot_dt_wt <- plot_dt[WT==T,.SD,,.SDcols = c("variable", "value")]

  # #Most abundant variant frequency
  # plot_dt_ma <- plot_dt[,.(value=max(value)),variable]

  #Fake single expected frequency: fake_single_counts = WT_count*error_rate/4
  plot_dt_fsingles <- data.table::copy(plot_dt_wt)
  plot_dt_fsingles[, value := ((value - 1)*error_rate/4 + 1)]
  plot_dt_fsingles[, Hamming_distance := factor(1)]

  #Fake double expected frequency distribution: fake_double_counts = real_single_counts*error_rate*2/4
  plot_dt_fdoubles <- plot_dt[get(fcols[1]) == 1]
  for(i in unique(plot_dt_fdoubles[,variable])){
    #Adjust observed single counts for outflow to doubles
    plot_dt_fdoubles[variable==i,real_value := value - 1 + (value - 1 - plot_dt_fsingles[variable==i,value - 1])*seq_length*error_rate]
    #Adjust observed single counts for inflow from WT
    plot_dt_fdoubles[variable==i,real_value := real_value - plot_dt_fsingles[variable==i,value - 1]]
    #Fake doubles
    plot_dt_fdoubles[variable==i,value := real_value*error_rate*2/4 + 1]
  }
  plot_dt_fdoubles <- plot_dt_fdoubles[,.(value = median(value)),variable]
  plot_dt_fdoubles[, Hamming_distance := factor(2)]

  #Remove real variants with hamming distance ==0
  plot_dt <- plot_dt[get(fcols[1])>0,]

  #Plot
  if(length(fcols)==1){
    #Plot without subplot stratification (Nham_nt)
    d <- ggplot2::ggplot(plot_dt, ggplot2::aes(value, color = Hamming_distance)) +
      ggplot2::geom_density() + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
      ggplot2::labs(x = "variant count + 1 (log scale)", y = "Density", title = title, color='Nucleotide\nHamming dist.') +
      # ggplot2::geom_vline(data = plot_dt_ma, ggplot2::aes(xintercept = value), linetype = 2) +
      ggplot2::geom_vline(data = plot_dt_wt, ggplot2::aes(xintercept = value), linetype = 2) +
      ggplot2::geom_vline(data = plot_dt_fsingles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
      ggplot2::geom_vline(data = plot_dt_fdoubles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
      ggplot2::facet_grid(Hamming_distance ~ variable, scales = "free")
    ggplot2::ggsave(output_file, width = width, height = height)
  }else{
    #Plot with subplot stratification (using second facet column, Nham_aa)
    plot_list <- list()
    plot_indices <- unlist(plot_dt[order(get(fcols[2]), decreasing = F),unique(.SD),,.SDcols = fcols[2]])
    for(i in plot_indices){
      this_plot_dt <- plot_dt[get(fcols[2])==i & get(fcols[1])<=max((i*3), 3) & get(fcols[1])>=i,]
      this_plot_dt[, Hamming_distance := factor(.SD[[1]], levels = 1:max((i*3), 3)),,.SDcols = fcols]
      plot_list[[as.character(i)]] <- ggplot2::ggplot(this_plot_dt, ggplot2::aes(value, color = Hamming_distance)) +
        ggplot2::geom_density() + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
        ggplot2::labs(x = "variant count + 1 (log scale)", y = "Density", title = paste0("Amino acid Hamming distance: ", i), color='Nucleotide\nHamming dist.') +
        # ggplot2::geom_vline(data = plot_dt_ma, ggplot2::aes(xintercept = value), linetype = 2) +
        ggplot2::geom_vline(data = plot_dt_wt, ggplot2::aes(xintercept = value), linetype = 2) +
        ggplot2::geom_vline(data = plot_dt_fsingles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
        ggplot2::geom_vline(data = plot_dt_fdoubles, ggplot2::aes(xintercept = value, color = Hamming_distance), linetype = 2) +
        ggplot2::facet_grid(Hamming_distance ~ variable, scales = "free")
    }
    #Combine plots
    temp_grid.arrange <- gridExtra::grid.arrange
    p <- do.call("temp_grid.arrange", c(plot_list, nrow = length(plot_indices)))
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
  }
}

