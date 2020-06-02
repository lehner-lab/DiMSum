
#' dimsum__error_model_qqplot
#'
#' Perform leave one out cross validation on replicates and generate QQ plot
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param norm_dt data.table of normalisation parameters (required)
#' @param error_dt data.table of error model parameters (required)
#' @param report_outpath fitness report output path (default:NULL)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__error_model_qqplot <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  norm_dt,
  error_dt,
  report_outpath = NULL
  ){

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")
  
  #Calculate error and fitness on training replicates and compare to fitness in leftout replicate
  for (i in seq_along(all_reps)) {
    
    training_reps <- paste0(strsplit(all_reps_str,"")[[1]][-i],collapse="")
    training_reps_num <- as.numeric(strsplit(training_reps,"")[[1]])
    NTreps <- nchar(training_reps)
    test_rep <- strsplit(all_reps_str,"")[[1]][i]
    training_data <- data.table::copy(input_dt)

    #Use error model parameters to calculate replicate-specific errors per variant
    for (j in seq_along(training_reps_num)) {
      if(dimsum_meta[["fitnessNormalise"]]){
        Corr <- matrix(unlist(norm_dt[,.SD,.SDcols = paste0("scale_", training_reps_num[j])]), ncol = 1, nrow = nrow(input_dt))
      }else{
        Corr <- matrix(1, ncol = 1, nrow = nrow(input_dt))
      }
      
      training_data[,paste0("error", training_reps_num[j]) := sqrt(Corr * rowSums(matrix(unlist(error_dt[parameter %in% c("input", "output") & rep == training_reps_num[j], mean_value]), nrow = .N, ncol = 2, byrow = T)/.SD) + 
        matrix(error_dt[parameter %in% c("reperror") & rep == training_reps_num[j], mean_value], nrow = .N, ncol = 1, byrow = T)),,
        .SDcols = c(
          grep(paste0("count_e", training_reps_num[j], "_s0"), names(training_data)), 
          grep(paste0("count_e", training_reps_num[j], "_s1"), names(training_data)))]
    }

    #Merged fitness values
    training_data[,fitness := rowSums(.SD[,1:NTreps]/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T) / 
                    rowSums(1/(.SD[,(NTreps+1):(2*NTreps)]^2),na.rm=T),
                  ,.SDcols = c(grep(paste0("^fitness[", training_reps, "]$"),names(training_data)),
                               grep(paste0("^error[", training_reps, "]$"),names(training_data)))]
    #Merged error
    training_data[,error := sqrt(1/rowSums(1/.SD^2)),,
                  .SDcols = grep(paste0("^error[", training_reps, "]$"),names(training_data))]
    
    #Predict test replicate ###
    #Only use variants that have non-NA fitness and error values
    training_data[,test_rep_ok := !is.na(.SD) & is.finite(unlist(.SD)) & all_reads == T & input_above_threshold == T,,.SDcols = paste0("fitness",test_rep)]

    Mi_full <- error_dt[parameter == "input" & rep == test_rep,mean(mean_value)]
    Mo_full <- error_dt[parameter == "output" & rep == test_rep,mean(mean_value)]
    A_full <- error_dt[parameter == "reperror" & rep == test_rep, mean_value]

    #Scaling factor from fitness normalization
    if(dimsum_meta[["fitnessNormalise"]]){
      fitness_scale <- norm_dt[,unlist(.SD),.SDcols = paste0("scale_",test_rep)]
    }else{
      fitness_scale <- 1
    }
    
    dt_zscore <- training_data[test_rep_ok==T,
        .(test_rep, Nham_nt,
            MioA_full = (fitness-unlist(.SD[,1]))/sqrt(error^2 + fitness_scale*(Mi_full/unlist(.SD[,3]) + Mo_full/unlist(.SD[,4])) + A_full)),
       ,.SDcols = c(paste0("fitness",test_rep),
                    paste0("cbe",test_rep),
                    names(training_data)[grep(paste0("count_e", test_rep, "_s0"), names(training_data))],
                    names(training_data)[grep(paste0("count_e", test_rep, "_s1"), names(training_data))])]

    if (i == 1) {
      dt_zscore_all = dt_zscore
    } else {
      dt_zscore_all = rbind(dt_zscore_all, dt_zscore)
    }
  }

  #Plot results for dataset
  dt_zscore_melt <- data.table::melt(dt_zscore_all[sample(.N, min(1e4, .N))], id.vars = c("test_rep", "Nham_nt"))
  #Repeat data for all replicates combined
  dt_zscore_melt_all <- data.table::copy(dt_zscore_melt)
  dt_zscore_melt_all[, test_rep := "All"]
  dt_zscore_melt <- rbind(dt_zscore_melt, dt_zscore_melt_all)
  #Plot colours
  plot_cols <- c('darkgrey', dimsum__gg_color_hue(length(all_reps)))
  names(plot_cols) <- c("All", all_reps)
  zrange = 3.5
  p <- ggplot2::ggplot(dt_zscore_melt, ggplot2::aes(sample = value, color = test_rep)) +
    ggplot2::geom_abline(lty = 2) + 
    ggplot2::geom_qq(geom = "line") +
    ggplot2::scale_colour_manual(name=c("Replicate"), values=plot_cols, guide='legend') +
    ggplot2::coord_cartesian(xlim = c(-zrange, zrange), ylim = c(-zrange, zrange)) +
    ggplot2::labs(x = 'Expected (Normal) theoretical quantiles', y = 'Observed data quantiles') +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(report_outpath,"dimsum_stage_fitness_report_1_errormodel_leaveoneout_qqplot.pdf"), p, width = 6, height = 5)
  ggplot2::ggsave(file.path(report_outpath,"dimsum_stage_fitness_report_1_errormodel_leaveoneout_qqplot.png"), p, width = 6, height = 5)
}
