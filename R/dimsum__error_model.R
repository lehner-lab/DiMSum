
#' dimsum__error_model
#'
#' Error model analysis including determination of count-based, replicate and over-sequencing error terms.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#'
#' @return error model data.table
#' @export
#' @import data.table
dimsum__error_model <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  report = TRUE,
  report_outpath = NULL
  ){

  message("Fit error model...")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")
  
  work_data <- input_dt[,.SD,,.SDcols = c("Nmut_nt","WT",
    grep(names(input_dt),pattern=paste0("e[", all_reps_str, "]_s0_b"),value=T),
    grep(names(input_dt),pattern=paste0("e[", all_reps_str, "]_s1_b"),value=T))]

  #Add up counts for biological output reps
  for(E in all_reps){
    idx <- names(work_data)[grep(names(work_data),pattern = paste0("e",E,"_s1_b"))]
    work_data[,paste0("count_e",E,"_s1") := rowSums(.SD),,.SDcols = idx]
    names(work_data)[grep(names(work_data),pattern = paste0("e",E,"_s0_b"))] <- paste0("count_e",E,"_s0")
  }
  work_data <- work_data[,.SD,.SDcols = c("Nmut_nt","WT",names(work_data)[grep(names(work_data),pattern="^count")])]

  ### Calculate fitness
  ###########################

  for(j in all_reps){
    wt_corr <- as.numeric(work_data[WT==T, log(.SD[,2]/.SD[,1]),,
      .SDcols = c(grep(paste0("count_e", j, "_s0"), names(work_data)), grep(paste0("count_e", j, "_s1"), names(work_data)))])
    work_data[, paste0("fitness",j) := log(.SD[,2]/.SD[,1]) - wt_corr,,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)), 
        grep(paste0("count_e", j, "_s1"), names(work_data)))]
    #Set infinite or undefined fitness values to NA
    work_data[is.nan(get(paste0("fitness",j))) | is.infinite(get(paste0("fitness",j))), paste0("fitness",j) := NA]
  }  

  #Flag variants that don't have reads in all input/output replicates
  work_data[, all_reads := rowSums(.SD > 0) == (2*nchar(all_reps_str)),,.SDcols = grep(paste0("count_e[", all_reps_str, "]_s[01]"), names(work_data))]

  ### Find input read threshold for full fitness range
  ###########################

  input_count_threshold <- work_data[all_reads == T,exp(-quantile(.SD,probs = 0.01,na.rm=T)),,.SDcols = grep("fitness",names(work_data))]
  #Define variants above threshold for later use
  work_data[,input_above_threshold := rowSums(.SD > input_count_threshold) == nchar(all_reps_str),,.SDcols = grep(paste0("count_e[", all_reps_str, "]_s0"),names(work_data))]
  
  #Plot input counts versus fitness with threshold
  if(report){
    X <- data.table::data.table()
    for(j in all_reps){
      if(length(X) == 0) {
        X = work_data[,.SD,,.SDcols = c(paste0("count_e", j, "_s0"),paste0("fitness",j))] 
        names(X) = c("input","fitness")
        X[,rep:=j]
      }else{
        Y = work_data[,.SD,,.SDcols = c(paste0("count_e", j, "_s0"),paste0("fitness",j))] 
        names(Y) = c("input","fitness")
        Y[,rep:=j]
        X = rbind(X,Y)
      }
    }
    d <- ggplot2::ggplot(X[input>0 & !is.na(fitness),],ggplot2::aes(input,fitness)) +
      ggplot2::geom_hex(size = 0.2, color = "lightgrey") +
      # ggplot2::scale_fill_viridis_c(trans="log10") +
      ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
      ggplot2::scale_x_log10() +
      ggplot2::geom_vline(xintercept = input_count_threshold, lty = 2) +
      ggplot2::geom_hline(yintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::facet_wrap(rep ~ .) + ggplot2::theme_bw() +
      ggplot2::labs(x = "Input variant count (log scale)", y = "Fitness", title = "Replicate fitness versus input variant counts")
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.pdf"), d, width = 6, height = 6)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.png"), d, width = 6, height = 6)
  }

  ### QC plots for fitness between replicates (before normalisation)
  ###########################

  if(report){
    #Fitness distributions (all replicates)
    X <- reshape2::melt(work_data[input_above_threshold == T & all_reads == T,.SD,.SDcols=grep("fitness",names(work_data))],measure.vars=grep("fitness",names(work_data),value=T))
    X[,replicate := gsub("fitness", "", as.character(variable)),variable]
    d <- ggplot2::ggplot(X,ggplot2::aes(value,color=replicate)) +
      ggplot2::geom_density() +
      ggplot2::labs(x = "Fitness", y = "Density", color = "replicate", title = "Before inter-replicate normalisation") + 
      ggplot2::geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::theme_bw()
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.pdf"), d, width = 6, height = 6)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.png"), d, width = 6, height = 6)

    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps)], 
      output_file = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.pdf"),
      xlab = "Fitness",
      ylab = "Fitness",
      title = "Before inter-replicate normalisation")
    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps)], 
      output_file = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.png"),
      xlab = "Fitness",
      ylab = "Fitness",
      title = "Before inter-replicate normalisation")
  }

  ### Calculate replicate normalisation parameters (scale and center/shift) to minimise inter-replicate differences
  ###########################

  #Fitness data
  F_data <- work_data[input_above_threshold == T & all_reads == T, as.matrix(.SD),.SDcols = grep(paste0("fitness[", all_reps_str, "]$"),names(work_data))]
  
  #Calculate replicate normalisation parameters using non-linear minimization
  set.seed(1603)
  p <- nlm(f = dimsum__replicate_fitness_deviation, p = rep(c(1,0), each = nchar(all_reps_str)), fitness_mat = F_data, all_reps = all_reps)[["estimate"]]
  
  #Normalise to first replicate (set scaling factor of first replicate to unity)
  p[1:nchar(all_reps_str)] <- p[1:nchar(all_reps_str)]/p[1]

  #Save normalisation model
  fitness_norm_model <- data.table(t(p))
  names(fitness_norm_model) <- c(paste0("scale_", all_reps), paste0("shift_", all_reps))
  write.table(round(fitness_norm_model, digits = 4),
    file = file.path(dimsum_meta[["tmp_path"]], "normalisationmodel.txt"), row.names = F, quote = F)
  
  #Wild-type correction such that mean(wild-type) = 0
  wt_corr <- work_data[WT == T, rowMeans((.SD + 
      unlist(fitness_norm_model[,.SD,,.SDcols = grep(paste0("shift_[", all_reps_str, "]$"), names(fitness_norm_model))])) * 
      unlist(fitness_norm_model[,.SD,,.SDcols = grep(paste0("scale_[", all_reps_str, "]$"), names(fitness_norm_model))])),,
    .SDcols = grep(paste0("fitness[", all_reps_str, "]$"), names(work_data))]
  
  #Normalize fitness
  for (j in all_reps){
    work_data[all_reads == T, paste0("fitness", j, "_norm") := (.SD + 
        unlist(fitness_norm_model[,.SD,,.SDcols = paste0("shift_", j)])) * 
        unlist(fitness_norm_model[,.SD,,.SDcols = paste0("scale_", j)]) - wt_corr,,
      .SDcols = paste0("fitness", j)]
  }

  if(report){
    #Fitness distributions (all replicates)
    X <- reshape2::melt(work_data[input_above_threshold == T & all_reads==T,.SD,.SDcols=grep("_norm",names(work_data))],measure.vars=grep("_norm",names(work_data),value=T))
    X[,replicate := gsub("_norm", "", gsub("fitness", "", as.character(variable))),variable]
    d <- ggplot2::ggplot(X,ggplot2::aes(value,color=replicate)) +
      ggplot2::geom_density() +
      ggplot2::labs(x = "Fitness", y = "Density", color = "replicate", title = "After inter-replicate normalisation") + 
      ggplot2::geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::theme_bw()
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.pdf"), d, width = 6, height = 6)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.png"), d, width = 6, height = 6)

    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps, "_norm")], 
      output_file = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.pdf"),
      xlab = "Fitness",
      ylab = "Fitness",
      title = "After inter-replicate normalisation")
    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps, "_norm")], 
      output_file = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.png"),
      xlab = "Fitness",
      ylab = "Fitness",
      title = "After inter-replicate normalisation")

    #Fitness replicate deviations (all replicates)
    X <- work_data[all_reads == T & input_above_threshold == T & Nmut_nt > 0, cbind(.SD - rowMeans(.SD), M = rowMeans(.SD)),
      .SDcols = grep(paste0("fitness[", all_reps_str, "]$"), names(work_data))]
    Xnorm <- work_data[all_reads == T & input_above_threshold == T & Nmut_nt > 0, cbind(.SD - rowMeans(.SD), M = rowMeans(.SD)),
      .SDcols = grep(paste0("fitness[", all_reps_str, "]_norm$"), names(work_data))]
    Y <- rbind(reshape2::melt(X, id.vars = "M"), reshape2::melt(Xnorm, id.vars = "M"))
    Y[, replicate := strsplit(as.character(variable), "_")[[1]][1], variable]
    Y[, normalised := grepl('norm', variable), variable]
    d <- ggplot2::ggplot(Y,ggplot2::aes(M, value, color = normalised)) +
      ggplot2::geom_smooth() +
      ggplot2::geom_hline(yintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::labs(x = "Mean fitness", y = "Deviation from mean fitness", title = "Fitness replicate deviations (all replicates)") + 
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(. ~ replicate)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicate_deviation_scatter.pdf"), d, width = 6, height = 6)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicate_deviation_scatter.png"), d, width = 6, height = 6)  
  }

  #Replace fitness with normalised fitness
  if(dimsum_meta[["fitnessNormalise"]]){
    for(j in all_reps){
      work_data[all_reads ==T,paste0("fitness", j) := .SD,
                ,.SDcols = paste0("fitness", j, "_norm")]
      # work_data[all_reads == T, paste0("fitness", j, "_norm") := NULL]
    }
  }
  
  ### Error model
  ###########################
  
  #First calculate count-based error per replicate
  for (j in all_reps) {
    wt_corr = as.numeric(work_data[WT==T,1/.SD[,2] + 1/.SD[,1],,
      .SDcols = c(grep(paste0("count_e", j, "_s0"),names(work_data)),grep(paste0("count_e", j, "_s1"),names(work_data)))])

    if(dimsum_meta[["fitnessNormalise"]]){
      work_data[,paste0("cbe",j) := sqrt(unlist(fitness_norm_model[,.SD,.SDcols = paste0("scale_", j)])) * sqrt(1/.SD[,2] + 1/.SD[,1] + wt_corr),,
        .SDcols = c(grep(paste0("count_e", j, "_s0"),names(work_data)),grep(paste0("count_e", j, "_s1"),names(work_data)))]
    }else{
      work_data[,paste0("cbe",j) := sqrt(1/.SD[,2] + 1/.SD[,1] + wt_corr),,
        .SDcols = c(grep(paste0("count_e", j, "_s0"),names(work_data)),grep(paste0("count_e", j, "_s1"),names(work_data)))]
    }
    #Set error of NA fitness values to NA
    work_data[is.na(get(paste0("fitness",j))),paste0("cbe",j) := NA]
  }

  #Then calculate density of data along mean count based error (for density dependent weighting)
  bins <- 50
  work_data[,mean_cbe := rowMeans(.SD),.SDcols = grep(paste0("^cbe[",all_reps_str,"]$"), names(work_data))]
  error_range <- seq(work_data[input_above_threshold == T & all_reads == T, log10(quantile(mean_cbe^2, probs=0.001))], 0, length.out = bins)
  D <- density(x = work_data[input_above_threshold == T & all_reads == T, log10(mean_cbe^2)],
    from = work_data[,log10(quantile(mean_cbe^2, probs = 0.001, na.rm = T))], to = 0, n = bins)
  work_data[,bin_error := findInterval(mean_cbe^2,vec = 10^error_range)]
  work_data[,bin_error_density := D$y[bin_error],bin_error]
  work_data[bin_error == 0, bin_error_density := work_data[bin_error > 0][which.min(bin_error), unique(bin_error_density)]]

  #Fit and write error model to file
  error_model <- dimsum__fit_error_model(
    dimsum_meta = dimsum_meta,
    input_dt = data.table::copy(work_data),
    all_reps = all_reps,
    norm_dt = fitness_norm_model)
  write.table(error_model, file = file.path(dimsum_meta[["tmp_path"]], "errormodel.txt"), row.names = F, quote = F)
  
  #Use error model parameters to calculate replicate-specific errors per variant
  for(j in all_reps){
    if(dimsum_meta[["fitnessNormalise"]]){
      Corr <- matrix(unlist(fitness_norm_model[,.SD,.SDcols = paste0("scale_", j)]), ncol = 1, nrow = nrow(work_data))
    }else{
      Corr <- matrix(1, ncol = 1, nrow = nrow(work_data))
    }

    work_data[,paste0("error", j) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("input", "output") & rep == j, mean_value]), nrow = .N, ncol = 2, byrow = T)/.SD) + 
      matrix(error_model[parameter %in% c("reperror") & rep == j, mean_value], nrow = .N, ncol = 1, byrow = T)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)), 
        grep(paste0("count_e", j, "_s1"), names(work_data)))]
  }

  ### Plot model parameters and fit
  ###########################

  if(report){
    #Add columns for upper/lower bounds of parameter estimates (for plotting)
    plot_error_model <- data.table::copy(error_model)
    plot_error_model[,upper := mean_value + sd_value]
    plot_error_model[,lower := mean_value - sd_value]
    plot_error_model[lower < 0,lower := mean_value]
    # print(plot_error_model)

    #For plot, calculate average variance of fitness values per bin
    bs_data <- work_data[Nmut_nt > 0 & input_above_threshold == T & all_reads ==T]
    bs_data[,v := apply(.SD,1,var),.SDcols = grep(paste0("^fitness[", all_reps_str, "]$"), names(bs_data))]
    bs_data[,bin_mean_var := mean(v), bin_error]
    bs_data[,bin_N := .N, bin_error]
    bs_data[,bin_mean_error := mean(mean_cbe^2), bin_error]
    #Average merged error per bin
    bs_data[,bin_pred_error := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error", all_reps)]
    melt_bs_data <- unique(data.table::melt(bs_data[,.SD,.SDcols = c("bin_pred_error","bin_error","bin_mean_error")], id.vars = c("bin_error","bin_mean_error")))
    
    #Plot1: input and output over-sequencing factor parameters +- sd
    a <- ggplot2::ggplot(plot_error_model[parameter %in% c("input","output")],
      ggplot2::aes(x=interaction(parameter, rep), mean_value, ymin = lower, ymax = upper, color = factor(rep), lty = parameter, shape = parameter)) +
      ggplot2::geom_pointrange(show.legend = F) +
      # ggplot2::scale_y_log10(limits = c(min(c(1, plot_error_model[parameter %in% c("input","output"), mean_value], plot_error_model[parameter %in% c("input", "output"), lower])),
      #   max(c(2.5, plot_error_model[parameter %in% c("input","output"), mean_value], plot_error_model[parameter %in% c("input", "output"), upper])))) +
      ggplot2::scale_y_log10() + ggplot2::theme_bw() + 
      ggplot2::labs(y = "Over-sequencing factor", x = "Replicate (input or output)")

    #Plot2: replicate error parameters +- sd
    b <- ggplot2::ggplot(plot_error_model[parameter == "reperror"], ggplot2::aes(x = factor(rep), y = sqrt(mean_value), ymin = sqrt(lower), ymax = sqrt(upper), color = factor(rep))) +
      ggplot2::geom_pointrange(show.legend = F) +
      # ggplot2::scale_y_log10(limits = c(min(plot_error_model[parameter == "reperror", sqrt(lower)]),
      #   max(c(0.1, plot_error_model[parameter == "reperror", sqrt(upper)])))) +
      ggplot2::scale_y_log10() + ggplot2::theme_bw() + 
      ggplot2::labs(y = "Replicate error", x = "Replicate")

    #Plot3: count based error against variance of fitness + fit
    temp_mean_error <- unique(bs_data[,.(bin_mean_error, bin_mean_var, bin_N)])
    temp_mean_error[, ymin := bin_mean_var * (1-2/bin_N)]
    temp_mean_error[, ymax := bin_mean_var * (1+2/bin_N)]
    temp_mean_error[ymin < bs_data[,quantile(v, 0.001)], ymin := bs_data[,quantile(v, 0.001)]]
    temp_mean_error[ymax > bs_data[,quantile(v, 0.999)], ymax := bs_data[,quantile(v, 0.999)]]
    c <- ggplot2::ggplot(bs_data[v>=bs_data[,quantile(v, 0.001)] & v<=bs_data[,quantile(v, 0.999)],]) +
      ggplot2::stat_binhex(ggplot2::aes(mean_cbe^2, v), bins = 100, size = 0.2, color = "lightgrey") +
      ggplot2::geom_pointrange(inherit.aes = F, data = temp_mean_error,
        ggplot2::aes(x = bin_mean_error, y = bin_mean_var, ymin = ymin, ymax = ymax), color = "orange") +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data,
        ggplot2::aes(bin_mean_error, value), lty = 2, size = 1, color = "black") +
      ggplot2::scale_y_log10(limits = c(bs_data[,quantile(v, 0.001)], bs_data[,quantile(v, 0.999)])) +
      ggplot2::scale_x_log10() +
      ggplot2::geom_abline(color = "darkgrey",lty=2) +
      ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
      ggplot2::coord_cartesian(xlim = c(min((bs_data[,mean_cbe])^2), 10^max(error_range))) +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = "Average count-based (Poisson) fitness error", y = "Variance of fitness values between replicates")

    #Combine plots
    p <- gridExtra::grid.arrange(a, b, c, nrow = 2,
      layout_matrix = rbind(c(1, 2), c(3, 3)), heights = c(1, 2), widths = c(2, 1))
    #Save
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_repspec.pdf"), p, width = 9, height = 9)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_repspec.png"), p, width = 9, height = 9)
  }

  message("Done")

  return(list(
    "error_model" = error_model,
    "norm_model" = fitness_norm_model))
}
