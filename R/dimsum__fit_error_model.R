
#' dimsum__fit_error_model
#'
#' Fit error model.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param norm_dt data.table of normalisation parameters (default:NULL)
#' @param Nbootstraps number of bootstraps (default:100)
#' @param maxN maximum sample size for each bootstrap (default:5000)
#' @param max_tries_per_fit maximum number of fitting attempts (default:20)
#' @param lower_rep replicate error lower bound (default:1e-4)
#'
#' @return Nothing
#' @export
dimsum__fit_error_model <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  norm_dt,
  Nbootstraps = 100,
  maxN = 10000,
  max_tries_per_fit = 20,
  lower_rep = 1e-4
  ){

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Fitness correction for normalisation (scaling) factors
  Fcorr <- NULL
  if(dimsum_meta[["fitnessNormalise"]]){
    Fcorr <- unlist(norm_dt[,.SD,.SDcols = grep(paste0("scale_(", all_reps_str, ")$"), names(norm_dt))])
  }

  #Set up fitting parameters
  Nreps <- length(all_reps)
  #Calculate all combinations of replicates of length >= 2
  idx <- list()
  for(i in (length(all_reps)):2){
    idx <- c(idx, combn(length(all_reps), i, function(x) list(x)))
  }

  #Subset to maximum 500 combinations (preferably highest order combinations)
  if(length(idx)>500){
    idx <- idx[1:500]
  }

  #Setup cluster
  clust <- parallel::makeCluster(dimsum_meta[['numCores']])
  #Set cluster seed
  parallel::clusterSetRNGStream(cl = clust, 1234567)
  #Fit
  parameters <- t(parallel::parSapply(
    cl = clust, 
    X = 1:Nbootstraps, 
    FUN = dimsum__fit_error_model_bootstrap, 
    input_dt = input_dt,
    all_reps = all_reps,
    idx = idx,
    maxN = maxN,
    max_tries_per_fit = max_tries_per_fit,
    lower_rep = lower_rep,
    Fcorr = Fcorr))
  #Stop cluster
  parallel::stopCluster(clust)
  
  #Return table with error model parameters
  error_model <- data.table(
    parameter = rep(c("input", "output", "reperror"), each = Nreps),
    rep = rep(all_reps, 3),
    mean_value = colMeans(parameters, na.rm = T),
    CI90_lower = apply(parameters, 2, quantile, na.rm = T, probs = 0.1),
    CI90_upper = apply(parameters, 2, quantile, na.rm = T, probs = 0.9),
    ensemble = sum(!is.na(parameters[,1])))
  return(error_model)
}
