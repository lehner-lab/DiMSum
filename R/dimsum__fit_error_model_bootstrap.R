
#' dimsum__fit_error_model_bootstrap
#'
#' Fit error model on random subsample of the input data.
#'
#' @param m bootstrap number; dummy variable (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param idx vector of replicate combinations (required)
#' @param maxN maximum size of subsample (required)
#' @param max_tries_per_fit maximum number of fitting attempts (required)
#' @param lower_rep replicate error lower bound (required)
#' @param Fcorr vector of fitness normalization factors (default:NULL)
#'
#' @return Nothing
#' @export
dimsum__fit_error_model_bootstrap <- function(
  m,
  input_dt,
  all_reps,
  idx,
  maxN,
  max_tries_per_fit,
  lower_rep,
  Fcorr = NULL
  ){
  require(data.table)

  Nreps <- length(all_reps)
  y <- NA
  counter <- 0
  while(!is.list(y) & counter < max_tries_per_fit){ #this is a routine to catch failed fitting attempts
    tryCatch({

      #Create data concatenation for all combinations of replicates (can be improved?!)
      bs_data <- input_dt[input_above_threshold == T & all_reads == T & Nham_nt > 0][sample(.N,min(.N,maxN),replace = T)]
      F_data_list <- list() #fitness data
      E_data_list <- list() #count based error data (for weighting of datapoints)
      C_data_list <- list() #inverse count data
      NR <- list() # simple counter for how many replicates are in each data row
      for(i in 1:length(idx)){
        F_data_list[[i]] <- bs_data[,.SD,.SDcols = grep(paste0("^fitness[", paste0(all_reps[idx[[i]]], collapse = ""), "]$"), names(input_dt))]
        E_data_list[[i]] <- bs_data[,.SD,.SDcols = grep(paste0("^cbe[", paste0(all_reps[idx[[i]]], collapse = ""), "]$"), names(input_dt))]
        C_data_list[[i]] <- bs_data[,1/.SD,.SDcols = grep(paste0("count_e[", paste0(all_reps[idx[[i]]], collapse = ""), "]_s"), names(input_dt))]
        NR[[i]] <- data.table(rep(length(idx[[i]]),nrow(bs_data)))
      }

      #make matrices with #replicates = #columns
      F_data <- as.matrix(rbindlist(F_data_list, fill = T)) 
      #Calculate variance from fitness data
      V_data <- apply(F_data,1,var,na.rm = T)

      #count based error weighting of data points (needed because deviation is proportional to expcectation)
      E_data <- as.matrix(rbindlist(E_data_list, fill = T))
      Ew <- rowMeans(E_data, na.rm = T)^2

      C_data <- as.matrix(rbindlist(C_data_list, fill = T))
      #Avoid NA for fitting
      C_data[is.na(C_data)] <- 0
      #Correct count terms for fitness normalization factors
      if(!is.null(Fcorr)){
        C_data = C_data * matrix(rep(Fcorr, 2), nrow = nrow(C_data), ncol = Nreps*2, byrow = T)
      }

      #weighting according to #variants with same mutation
      Dw <- bs_data$error_model_weighting
      #how many replicates are in the replicate combinations
      NRT <- as.matrix(rbindlist(NR))
      
      #Binary variable to avoid replicate error fitting for replicates that are not present in a data row
      BV <- F_data
      BV[!is.na(BV)] <- 1
      BV[is.na(BV)] <- 0
      
      #fit formula
      fit_formula <- as.formula(paste0(
        "V_data ~ (",
        paste0("BV[,",
          1:Nreps,
          "] * (p[",
          1:Nreps,
          "] * C_data[,",
          1:Nreps,
          "] +  p[",
          Nreps + 1:Nreps,
          "] * C_data[,",
          Nreps + 1:Nreps,
          "] + p[",
          2*Nreps + 1:Nreps,
          "])",collapse = " + "),
        ") / NRT"))

      #Fitting
      y <- nls(formula = fit_formula,
            start = list(p = c(
              rep(1, each = Nreps*2) + 10^rnorm(2*Nreps), 
              rep(0, each = Nreps) + 0.01*10^rnorm(1*Nreps))),
            lower = c(rep(1, each = Nreps*2), rep(lower_rep, each = Nreps)),
            weights = 1 / Ew / Dw,
            algorithm = "port")
    }, error = function(cond){})
    counter <- counter + 1}
    
  #If fit was succesful after less than max_tries_per_fit attempts, transfer parameters to p
  if(counter < max_tries_per_fit){
    p <- summary(y)[["parameters"]][,1]
  }else{ #otherwise move on to next fit
    p <- rep(NA, 3*Nreps)
  }
  return(p)
}
