
#' dimsum__bayesian_double_fitness
#'
#' Calculating Bayesian double mutant fitness estimates.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param doubles_dt doubles data.table (required)
#' @param singles_dt singles data.table (required)
#' @param wt_dt WT data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param min_input_read_count_doubles minimum input read count for doubles used to derive prior for Bayesian doubles correction (required)
#' @param lam_d Poisson distribution for score likelihood (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__bayesian_double_fitness <- function(
  dimsum_meta,
  doubles_dt,
  singles_dt,
  wt_dt,
  all_reps,
  min_input_read_count_doubles,
  lam_d,
  report = TRUE,
  report_outpath = NULL
  ){

  message("Calculating Bayesian double mutant fitness estimates...")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #Plot fitness densities for different mean count bins (replicate 1)
  doubles_dt[,fitness_for_bins := .SD[[1]],,.SDcols = paste0("fitness",all_reps[1],"_uncorr")]
  if(report){
    d <- ggplot2::ggplot(doubles_dt[between(bin_count,2,8) & !is.infinite(fitness_for_bins)],ggplot2::aes(fitness_for_bins,..density..,color=factor(bin_count))) +
      ggplot2::geom_density(adjust=1)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework1.png"), d, width = 7, height = 5)
  }

  #Plot fitness densities for mean counts greater/less than 50 (replicate 1)
  if(report){
    d <- ggplot2::ggplot(doubles_dt[!is.infinite(fitness_for_bins) & !is.na(fitness_for_bins)],ggplot2::aes(fitness_for_bins,..density..,color=counts_for_bins >= min_input_read_count_doubles)) +
      ggplot2::geom_density(adjust=1)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework2.png"), d, width = 7, height = 5)
  }
  #>> try to estimate what fitness scores are for variants with low sequence coverage
  # use double mutants with variants >= min_input_read_count_doubles counts 

  ## calculate posterior double mutant fitness based on prior from single mutants
  postpois_conditioned_singleF <- function(i){  
    require(data.table)
    count_in <- double_data[i,count_in]
    count_out <- double_data[i,count_out]
    lam_in <- exp(seq(floor(log(count_in+0.1)-max(c(0.5,1/log10(count_in+1.75)))),(log(count_in+0.1)+max(c(0.5,1/log10(count_in+1.75)))),lam_d))
    lam_out <- exp(seq(floor(log(count_out+0.1)-max(c(0.5,1/log10(count_out+1.75)))),(log(count_out+0.1)+max(c(0.5,1/log10(count_out+1.75)))),lam_d))
    lam_low <- range(log(lam_out))[1] - range(log(lam_in))[2]
    lam_high <- range(log(lam_out))[2] - range(log(lam_in))[1]
    idx <- row(matrix(NA,nrow=length(lam_out),ncol=length(lam_in))) - col(matrix(NA,nrow=length(lam_out),ncol=length(lam_in)))
    likelihood <- sapply(split(outer(dpois(count_out,lambda = lam_out),dpois(count_in,lambda = lam_in)),idx),sum)
    score_prior <- density(score_prior_cond[,.(fdist = sqrt((double_data[i,F1]-F1)^2+(double_data[i,F2]-F2)^2),F)][
      order(fdist)][1:Nneighbours,F],
      from = (lam_low-wt_corr),
      to = (lam_high-wt_corr),
      n = as.integer(as.character(round((lam_high-lam_low)/lam_d + 1)))) #super weird bug
    posterior <- score_prior$y*likelihood
    
    moments <- list()
    moments[1] <- weighted.mean(x = score_prior$x,w = posterior)
    moments[2] <- sqrt(sum(( moments[[1]]-score_prior$x)^2 * posterior)/
                        sum(posterior))
    return(moments)
  }

  # Setup cluster
  clust <- parallel::makeCluster(dimsum_meta[['numCores']]) #This line will take time

  #Calculate conditional fitness and sigma
  for (E in all_reps) {
    #wildtype "correction" to calculate scores
    wt_corr <- wt_dt[,log(unlist(.SD[,1]) / unlist(.SD[,2])),,.SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    #data for prior calculation
    double_data <- doubles_dt[,.(Pos1,Mut1,Pos2,Mut2,count_in = unlist(.SD[,1]),count_out = unlist(.SD[,2]),
                             F = unlist(.SD[,3])),,
                          .SDcols = c(paste0("count_e",E,"_s0"),paste0("count_e",E,"_s1"),paste0("fitness",E,"_uncorr"))]
    # double_data = merge(double_data,singles_dt[,.(Pos,Mut,F1 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    double_data <- merge(double_data,singles_dt[!is.na(singles_dt[,paste0("fitness",E)]) & is.reads0==T,.(Pos,Mut,F1 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    # double_data = merge(double_data,singles_dt[,.(Pos,Mut,F2 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
    double_data <- merge(double_data,singles_dt[!is.na(singles_dt[,paste0("fitness",E)]) & is.reads0==T,.(Pos,Mut,F2 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
    
    Nneighbours <- 500
    score_prior_cond <- double_data[count_in >= min_input_read_count_doubles & F > -Inf & F1 > -Inf & F2 > -Inf]

    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("double_data","lam_d","wt_corr","score_prior_cond","Nneighbours"), envir = environment())

    #posterior fitness conditioned on single fitness
    t<-proc.time()
    helper <- parallel::parSapply(clust,X = 1:nrow(double_data), postpois_conditioned_singleF)
    print(proc.time()-t)
    helper1 <- matrix(unlist(helper),nrow=2)
    double_data[,paste0("fitness",E,"_cond") := helper1[1,]]
    double_data[,paste0("sigma",E,"_cond") := helper1[2,]]
    doubles_dt <- merge(doubles_dt, double_data[,.SD,,.SDcols = c("Pos1", "Pos2", "Mut1", "Mut2", paste0("fitness",E,"_cond"), paste0("sigma",E,"_cond"))], by = c("Pos1", "Pos2", "Mut1", "Mut2"), all.x = T)
  }
  parallel::stopCluster(clust)

  if(report){
    #Scatterplot matrix - singles
    if(dimsum_meta[["sequenceType"]]=="coding"){
      d <- GGally::ggpairs(singles_dt[Nmut_aa==1 & apply(abs(singles_dt[,.SD,,.SDcols = paste0("fitness",all_reps)])==Inf, 1, sum)==0,
        grep(names(singles_dt),pattern="fitness"),with=F],
        upper=list(continuous = "cor"))
    }else{
      d <- GGally::ggpairs(singles_dt[apply(abs(singles_dt[,.SD,,.SDcols = paste0("fitness",all_reps)])==Inf, 1, sum)==0,
        grep(names(singles_dt),pattern="fitness"),with=F],
        upper=list(continuous = "cor"))
    }
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework_scattermatrix_singles.png"), d, width = 10, height = 10)

    #Scatterplot matrix - doubles, uncorrected
    set.seed(1)
    d <- GGally::ggpairs(doubles_dt[apply(abs(doubles_dt[,.SD,,.SDcols = paste0("fitness",all_reps,"_uncorr")])==Inf, 1, sum)==0
      ][sample(x = .N,min(c(.N,1000))),grep(names(doubles_dt),pattern=paste0("fitness[", all_reps_str, "]_uncorr")),with=F],
      upper=list(continuous = "cor"))
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework_scattermatrix_doubles_uncorr.png"), d, width = 10, height = 10)

    #Scatterplot matrix - doubles, conditional
    set.seed(1)
    d <- GGally::ggpairs(doubles_dt[apply(abs(doubles_dt[,.SD,,.SDcols = paste0("fitness",all_reps,"_cond")])==Inf, 1, sum)==0
      ][sample(x = .N,min(c(.N,1000))),grep(names(doubles_dt),pattern=paste0("fitness[", all_reps_str, "]_cond")),with=F],
      upper=list(continuous = "cor"))
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_bayesian_framework_scattermatrix_doubles_cond.png"), d, width = 10, height = 10)
  }
  
  message("Done")

  return(doubles_dt)

}
