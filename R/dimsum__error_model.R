
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

  #Return NULL model if no fit required
  if(!dimsum_meta[["fitnessErrorModel"]]){
    return(list(
      "error_model" = NULL,
      "norm_model" = NULL))
  }

  dimsum__status_message("Fit error model...\n")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")
  
  work_data <- input_dt[,.SD,.SDcols = c("Nham_nt","WT",names(input_dt)[grep(names(input_dt),pattern="^count")])]

  ### Calculate fitness
  ###########################

  for(j in all_reps){
    work_data[, paste0("fitness",j) := log(.SD[,2]/.SD[,1]),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)), 
        grep(paste0("count_e", j, "_s1"), names(work_data)))]
    #Set infinite or undefined fitness values to NA
    work_data[is.nan(get(paste0("fitness",j))) | is.infinite(get(paste0("fitness",j))), paste0("fitness",j) := NA]
  }  

  #Flag variants that don't have reads in all input/output replicates
  work_data[, all_reads := rowSums(.SD > 0) == (2*nchar(all_reps_str)),,.SDcols = grep(paste0("count_e[", all_reps_str, "]_s[01]"), names(work_data))]

  #Check if WT data present in all input/output replicates
  if(work_data[WT == T & all_reads == T,.N]==0){
    input_dt[, all_reads := rowSums(.SD > 0) == (2*nchar(all_reps_str)),,.SDcols = grep(paste0("count_e[", all_reps_str, "]_s[01]"), names(input_dt))]
    input_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("count_e[", all_reps_str, "]_s0"), names(input_dt))]
    dimsum__status_message(paste0("WT variant has zero count in at least one input/output replicate. Did you mean to specify one of the following?\n"))
    if(dimsum_meta[["sequenceType"]]=="coding" & dimsum_meta[["mixedSubstitutions"]]){
      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(aa_seq, all_reads, mean_count)])
    }else if(dimsum_meta[["sequenceType"]]=="coding"){
      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), aa_seq, all_reads, mean_count)])
    }else{
      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
    }
    stop(paste0("Cannot proceed with error modelling: WT variant has zero count in at least one input/output replicate"), call. = FALSE)
  }

  ### Find input read threshold for full fitness range
  ###########################

  input_count_threshold <- work_data[all_reads == T,exp(-quantile(.SD,probs = 0.01,na.rm=T)),,.SDcols = grep("fitness",names(work_data))]
  #Define variants above threshold for later use
  work_data[,input_above_threshold := rowSums(.SD > input_count_threshold) == nchar(all_reps_str),,.SDcols = grep(paste0("count_e[", all_reps_str, "]_s0"),names(work_data))]

  #Correct for WT fitness
  for(j in all_reps){
    wt_corr <- as.numeric(work_data[WT==T, .SD,,.SDcols = paste0("fitness",j)])
    work_data[, paste0("fitness",j) := .SD - wt_corr,,.SDcols = paste0("fitness",j)]
  }  

  #Plot input counts versus fitness with threshold
  if(report){
    X <- data.table::data.table()
    for(j in all_reps){
      if(length(X) == 0) {
        X = work_data[,.SD,,.SDcols = c(paste0("count_e", j, "_s0"),paste0("fitness",j))] 
        names(X) = c("input","fitness")
        X[,rep := paste0("replicate", j)]
      }else{
        Y = work_data[,.SD,,.SDcols = c(paste0("count_e", j, "_s0"),paste0("fitness",j))] 
        names(Y) = c("input","fitness")
        Y[,rep := paste0("replicate", j)]
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
      ggplot2::labs(x = "Input variant count (log scale)", y = "Fitness")#, title = "Replicate fitness versus input variant counts")
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.pdf"), d, width = 6, height = 6)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.png"), d, width = 6, height = 6)
  }

  #Check if sufficient data remains above input read threshold (and data present in all input/output replicates)
  #Define minimum #variants for fitting = 10 x number of fitted parameters (3 x #reps)
  min_n_variants <- 10 * 3 * length(all_reps)
  if(work_data[input_above_threshold == T & all_reads == T & Nham_nt > 0,.N] < min_n_variants){
    stop(paste0("Cannot proceed with error modelling: insufficent number of variants satisfying full fitness range"), call. = FALSE)
  }

  ### QC plots for fitness between replicates (before normalisation)
  ###########################

  if(report){
    #Fitness distributions (all replicates)
    X <- reshape2::melt(work_data[input_above_threshold == T & all_reads == T,.SD,.SDcols=grep("fitness",names(work_data))],measure.vars=grep("fitness",names(work_data),value=T))
    X[,replicate := gsub("fitness", "", as.character(variable)),variable]
    d <- ggplot2::ggplot(X,ggplot2::aes(value,color=replicate)) +
      ggplot2::geom_density() +
      ggplot2::labs(x = "Fitness", y = "Density", color = "replicate") +#, title = "Before inter-replicate normalisation") + 
      ggplot2::geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::theme_bw()
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.pdf"), d, width = 6, height = 4)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.png"), d, width = 6, height = 4)

    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps)], 
      output_file = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.pdf"),
      xlab = "Fitness",
      ylab = "Fitness")#, title = "Before inter-replicate normalisation")
    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps)], 
      output_file = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.png"),
      xlab = "Fitness",
      ylab = "Fitness")#, title = "Before inter-replicate normalisation")
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
      ggplot2::labs(x = "Fitness", y = "Density", color = "replicate") +#, title = "After inter-replicate normalisation") + 
      ggplot2::geom_vline(xintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::theme_bw()
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.pdf"), d, width = 6, height = 4)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.png"), d, width = 6, height = 4)

    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps, "_norm")], 
      output_file = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.pdf"),
      xlab = "Fitness",
      ylab = "Fitness")#, title = "After inter-replicate normalisation")
    #Fitness correlations (all replicates)
    dimsum__ggpairs_binhex(
      input_dt = work_data[all_reads == T,.SD,,.SDcols = paste0("fitness", all_reps, "_norm")], 
      output_file = file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.png"),
      xlab = "Fitness",
      ylab = "Fitness")#, title = "After inter-replicate normalisation")

    #Fitness replicate deviations (all replicates)
    X <- work_data[all_reads == T & input_above_threshold == T & Nham_nt > 0, cbind(.SD - rowMeans(.SD), M = rowMeans(.SD)),
      .SDcols = grep(paste0("fitness[", all_reps_str, "]$"), names(work_data))]
    Xnorm <- work_data[all_reads == T & input_above_threshold == T & Nham_nt > 0, cbind(.SD - rowMeans(.SD), M = rowMeans(.SD)),
      .SDcols = grep(paste0("fitness[", all_reps_str, "]_norm$"), names(work_data))]
    Y <- rbind(reshape2::melt(X, id.vars = "M"), reshape2::melt(Xnorm, id.vars = "M"))
    Y[, replicate := strsplit(as.character(variable), "_")[[1]][1], variable]
    Y[, normalised := grepl('norm', variable), variable]
    d <- ggplot2::ggplot(Y,ggplot2::aes(M, value, color = normalised)) +
      ggplot2::geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
      ggplot2::geom_hline(yintercept = 0, lty = 2, color = "darkgrey") +
      ggplot2::labs(x = "Mean fitness", y = "Deviation from mean fitness") +#, title = "Fitness replicate deviations (all replicates)") + 
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

  #weight variants according to how many other variants with same # of mutations are around
  work_data[, error_model_weighting := sqrt(max(.N, sqrt(nrow(work_data)))), Nmut]

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

    #Full model
    work_data[,paste0("error", j) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("input", "output") & rep == j, mean_value]), nrow = .N, ncol = 2, byrow = T)/.SD) + 
      matrix(error_model[parameter %in% c("reperror") & rep == j, mean_value], nrow = .N, ncol = 1, byrow = T)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)), 
        grep(paste0("count_e", j, "_s1"), names(work_data)))]

    #Just replicate error
    work_data[,paste0("error_reponly", j) := sqrt(matrix(error_model[parameter %in% c("reperror") & rep == j, mean_value], nrow = .N, ncol = 1, byrow = T))]

    #Just input over-sequencing factor
    work_data[,paste0("error_inputosf", j) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("input") & rep == j, mean_value]), nrow = .N, ncol = 1, byrow = T)/.SD)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)))]

    #Just input (over-sequencing factor==1)
    work_data[,paste0("error_inputonly", j) := sqrt(Corr * rowSums(matrix(1, nrow = .N, ncol = 1, byrow = T)/.SD)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(work_data)))]

    #Just output over-sequencing factor
    work_data[,paste0("error_outputosf", j) := sqrt(Corr * rowSums(matrix(unlist(error_model[parameter %in% c("output") & rep == j, mean_value]), nrow = .N, ncol = 1, byrow = T)/.SD)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s1"), names(work_data)))]

    #Just output (over-sequencing factor==1)
    work_data[,paste0("error_outputonly", j) := sqrt(Corr * rowSums(matrix(1, nrow = .N, ncol = 1, byrow = T)/.SD)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s1"), names(work_data)))]

  }

  ### Plot model parameters and fit
  ###########################

  if(report){
    #Add columns for upper/lower bounds of parameter estimates (for plotting)
    plot_error_model <- data.table::copy(error_model)
    # print(plot_error_model)

    #Error model limits
    rep_error_intercept <- mean(plot_error_model[parameter=="reperror",mean_value])
    mult_error_slope <- median(plot_error_model[parameter!="reperror",mean_value])

    #For plot, calculate average variance of fitness values per bin
    bs_data <- work_data[Nham_nt > 0 & input_above_threshold == T & all_reads ==T]
    bs_data[,v := apply(.SD,1,var),.SDcols = grep(paste0("^fitness[", all_reps_str, "]$"), names(bs_data))]
    bs_data[,bin_mean_var := mean(v), bin_error]
    bs_data[,bin_N := .N, bin_error]
    bs_data[,bin_mean_error := mean(mean_cbe^2), bin_error]
    #Average merged error per bin
    bs_data[,bin_pred_error := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error", all_reps)]
    bs_data[,bin_pred_error_reponly := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_reponly", all_reps)]
    bs_data[,bin_pred_error_inputosf := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_inputosf", all_reps)]
    bs_data[,bin_pred_error_inputonly := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_inputonly", all_reps)]
    bs_data[,bin_pred_error_outputosf := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_outputosf", all_reps)]
    bs_data[,bin_pred_error_outputonly := mean(rowMeans(.SD^2, na.rm = T)), bin_error, .SDcols = paste0("error_outputonly", all_reps)]
    melt_bs_data <- unique(data.table::melt(bs_data[,.SD,.SDcols = c(
      "bin_pred_error",
      "bin_pred_error_reponly",
      "bin_pred_error_inputosf",
      "bin_pred_error_inputonly",
      "bin_pred_error_outputosf",
      "bin_pred_error_outputonly",
      "bin_error",
      "bin_mean_error")], id.vars = c("bin_error","bin_mean_error")))
    
    #Plot1: input and output over-sequencing factor parameters +- sd
    a <- ggplot2::ggplot(plot_error_model[parameter %in% c("input","output")],
      ggplot2::aes(x=interaction(parameter, rep), mean_value, ymin = CI90_lower, ymax = CI90_upper, color = factor(rep))) +
      ggplot2::geom_pointrange(show.legend = F) +
      # ggplot2::scale_y_log10(limits = c(min(c(1, plot_error_model[parameter %in% c("input","output"), mean_value], plot_error_model[parameter %in% c("input", "output"), CI90_lower])),
      #   max(c(2.5, plot_error_model[parameter %in% c("input","output"), mean_value], plot_error_model[parameter %in% c("input", "output"), CI90_upper])))) +
      ggplot2::scale_y_log10() + ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(y = "Multiplicative\nerror terms", x = "Replicate (input or output)")

    #Plot2: replicate error parameters +- sd
    b <- ggplot2::ggplot(plot_error_model[parameter == "reperror"], ggplot2::aes(x = factor(rep), y = sqrt(mean_value), ymin = sqrt(CI90_lower), ymax = sqrt(CI90_upper), color = factor(rep))) +
      ggplot2::geom_pointrange(show.legend = F) +
      # ggplot2::scale_y_log10(limits = c(min(plot_error_model[parameter == "reperror", sqrt(CI90_lower)]),
      #   max(c(0.1, plot_error_model[parameter == "reperror", sqrt(CI90_upper)])))) +
      ggplot2::scale_y_log10() + ggplot2::theme_bw() + 
      ggplot2::labs(y = "Additive\nerror terms", x = "Replicate")

    #Plot3: count based error against variance of fitness + fit
    plot_cols <- dimsum__gg_color_hue(6)
    temp_mean_error <- unique(bs_data[,.(bin_mean_error, bin_mean_var, bin_N)])
    temp_mean_error[, ymin := bin_mean_var * (1-2/bin_N)]
    temp_mean_error[, ymax := bin_mean_var * (1+2/bin_N)]
    temp_mean_error[ymin < bs_data[,quantile(v, 0.001)], ymin := bs_data[,quantile(v, 0.001)]]
    temp_mean_error[ymax > bs_data[,quantile(v, 0.999)], ymax := bs_data[,quantile(v, 0.999)]]
    c <- ggplot2::ggplot(bs_data[v>=bs_data[,quantile(v, 0.001)] & v<=bs_data[,quantile(v, 0.999)],]) +
      ggplot2::stat_binhex(ggplot2::aes(mean_cbe^2, v), bins = 100, size = 0.2, color = "lightgrey") +
      ggplot2::geom_abline(linetype = 2, size = 1) +
      # ggplot2::geom_abline(intercept = log10(mult_error_slope), linetype = 2, size = 1, color = plot_cols[2]) +
      # ggplot2::geom_hline(yintercept = rep_error_intercept, linetype = 2, size = 1, color = plot_cols[2]) +
      ggplot2::geom_pointrange(inherit.aes = F, data = temp_mean_error,
        ggplot2::aes(x = bin_mean_error, y = bin_mean_var, ymin = ymin, ymax = ymax), color = plot_cols[5]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error",],
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[1]) +
      ggplot2::scale_y_log10(limits = c(bs_data[,quantile(v, 0.001)], bs_data[,quantile(v, 0.999)])) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
      ggplot2::coord_cartesian(xlim = c(min((bs_data[,mean_cbe])^2), 10^max(error_range))) +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = "Count-based error estimate", y = "Variance of replicate fitness scores")

    #Plot4: error model decomposition
    d <- ggplot2::ggplot(melt_bs_data[variable=="bin_pred_error",]) +
      ggplot2::geom_line(inherit.aes = F,
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[1]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_reponly",],
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[3]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_inputosf",],
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[4]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_inputonly",],
        ggplot2::aes(bin_mean_error, value), size = 0.3, color = plot_cols[4], linetype = 2) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_outputosf",],
        ggplot2::aes(bin_mean_error, value), size = 1, color = plot_cols[6]) +
      ggplot2::geom_line(inherit.aes = F, data = melt_bs_data[variable=="bin_pred_error_outputonly",],
        ggplot2::aes(bin_mean_error, value), size = 0.3, color = plot_cols[6], linetype = 2) +
      ggplot2::scale_y_log10(limits = c(bs_data[,quantile(v, 0.001)], bs_data[,quantile(v, 0.999)])) +
      ggplot2::scale_x_log10() +
      ggplot2::coord_cartesian(xlim = c(min((bs_data[,mean_cbe])^2), 10^max(error_range))) +
      ggplot2::theme_bw() + 
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 8),
        axis.title.y = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(size = 7),
        axis.text.y = ggplot2::element_text(size = 7)) +
      ggplot2::labs(x = "Count-based error estimate", y = "Variance of replicate fitness scores")

    #Plot5: legend
    legend_names <- c(
        "Empirical variance (average per bin)", 
        "Poisson variance (null expectation)", 
        "Full model variance",
        "Input multiplicative terms only\n(and null expectation)",
        "Output multiplicative terms only\n(and null expectation)",
        "Additive terms only")
    legend_data <- data.frame(
      x = rep(1:10, times = 6),
      y = rep(1:10, times = 6), 
      variance_type = rep(legend_names, each = 10))
    legend_data[, "variance_type"] <- factor(legend_data[, "variance_type"], levels = legend_names)
    legend_cols <- c(
      plot_cols[5],
      "black",
      plot_cols[1],
      plot_cols[4],
      plot_cols[6],
      plot_cols[3])
    names(legend_cols) <- legend_names
    e <- ggplot2::ggplot(data = legend_data, ggplot2::aes(x,y)) + 
      ggplot2::geom_point(ggplot2::aes(color=variance_type)) +
      ggplot2::geom_line(ggplot2::aes(color=variance_type)) +
      ggplot2::scale_colour_manual(name=c("Legend"), values=legend_cols, guide='legend') +
      ggplot2::guides(
        colour = ggplot2::guide_legend(override.aes = list(
          linetype=c(NA,2,1,1,1,1), 
          size=c(2,0.75,0.75,0.75,0.75,0.75), 
          shape=c(19, NA, NA, NA, NA, NA)))) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = 9))
    e <- cowplot::get_legend(e)

    #Combine plots
    p <- gridExtra::grid.arrange(a, b, c, d, e, nrow = 3,
      layout_matrix = rbind(c(1, 2), c(3, 4), c(3, 5)), heights = c(1, 1, 1), widths = c(2, 0.8))
    #Save
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_repspec.pdf"), p, width = 9, height = 7)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_errormodel_repspec.png"), p, width = 9, height = 7)
  }

  dimsum__status_message("Done\n")

  #Render report
  dimsum__render_report(dimsum_meta = dimsum_meta)

  return(list(
    "error_model" = error_model,
    "norm_model" = fitness_norm_model))
}
