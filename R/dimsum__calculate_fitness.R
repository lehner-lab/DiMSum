
#' dimsum__calculate_fitness
#'
#' Calculate fitness and error.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param error_model_dt error model data.table (required)
#' @param norm_model_dt normalisation model data.table (required)
#' @param verbose whether or not to print status messages (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__calculate_fitness <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  error_model_dt,
  norm_model_dt,
  verbose = T
  ){

  if(verbose){message("Calculating fitness and error...")}

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #Calculate fitness
  for(j in all_reps){
    wt_corr <- as.numeric(input_dt[WT==T, log(.SD[,2]/.SD[,1]),,
      .SDcols = c(grep(paste0("count_e", j, "_s0"), names(input_dt)), grep(paste0("count_e", j, "_s1"), names(input_dt)))])
    input_dt[, paste0("fitness", j, "_uncorr") := log(.SD[,2]/.SD[,1]) - wt_corr,,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(input_dt)),
        grep(paste0("count_e", j, "_s1"), names(input_dt)))]
    #Set infinite or undefined fitness values to NA
    input_dt[is.nan(get(paste0("fitness",j, "_uncorr"))) | is.infinite(get(paste0("fitness",j, "_uncorr"))), paste0("fitness",j, "_uncorr") := NA]
  }  

  #Normalize fitness if necessary
  if(dimsum_meta[["fitnessNormalise"]]){
    #Wild-type correction such that mean(wild-type) = 0
    wt_corr <- input_dt[WT == T, rowMeans((.SD + 
        unlist(norm_model_dt[,.SD,,.SDcols = grep(paste0("shift_[", all_reps_str, "]$"), names(norm_model_dt))])) * 
        unlist(norm_model_dt[,.SD,,.SDcols = grep(paste0("scale_[", all_reps_str, "]$"), names(norm_model_dt))])),,
      .SDcols = grep(paste0("fitness[", all_reps_str, "]_uncorr$"), names(input_dt))]
    #Shift and scale
    for(j in all_reps){
      input_dt[, paste0("fitness", j, "_uncorr") := (.SD + 
          unlist(norm_model_dt[,.SD,,.SDcols = paste0("shift_", j)])) * 
          unlist(norm_model_dt[,.SD,,.SDcols = paste0("scale_", j)]) - wt_corr,,
        .SDcols = paste0("fitness", j, "_uncorr")]
    }
  }

  #Use error model parameters to calculate replicate-specific errors per variant
  for(j in all_reps){
    if(dimsum_meta[["fitnessNormalise"]]){
      Corr <- matrix(unlist(norm_model_dt[,.SD,.SDcols = paste0("scale_", j)]), ncol = 1, nrow = nrow(input_dt))
    }else{
      Corr <- matrix(1, ncol = 1, nrow = nrow(input_dt))
    }

    input_dt[, paste0("sigma", j, "_uncorr") := sqrt(Corr * rowSums(matrix(unlist(error_model_dt[parameter %in% c("input", "output") & rep == j, mean_value]), nrow = .N, ncol = 2, byrow = T)/.SD) + 
      matrix(error_model_dt[parameter %in% c("reperror") & rep == j, mean_value], nrow = .N, ncol = 1, byrow = T)),,
      .SDcols = c(
        grep(paste0("count_e", j, "_s0"), names(input_dt)), 
        grep(paste0("count_e", j, "_s1"), names(input_dt)))]
    #Set error of NA fitness values to NA
    input_dt[is.na(get(paste0("fitness", j, "_uncorr"))),paste0("sigma", j, "_uncorr") := NA]
  }

  #Remove unnecessary columns
  output_dt <- input_dt[,.SD,merge_seq,.SDcols = c(
    "aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","STOP","STOP_readthrough",names(input_dt)[grep(names(input_dt),pattern="^count|^fitness|^sigma")])]

  #Remove variants without fitness estimates in any replicates
  fitness_cols <- names(output_dt)[grep(paste0("fitness[", all_reps_str, "]_uncorr"), names(output_dt))]
  output_dt <- output_dt[!is.nan(rowMeans(output_dt[,fitness_cols,with=F], na.rm = T))]

  #Calculate mean input counts
  output_dt[,mean_count := rowMeans(.SD, na.rm = T),,.SDcols = paste0("count_e", all_reps, "_s0")]

  if(verbose){message("Done")}

  return(output_dt)

}
