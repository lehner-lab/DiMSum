
#' dimsum__add_dropout_pseudocount
#'
#' Add pseudocount to output samples with dropout (zero read counts).
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param verbose whether or not to print status messages (default: TRUE)
#'
#' @return data.table with pseudocounts added to output samples with dropout
#' @export
#' @import data.table
dimsum__add_dropout_pseudocount <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  verbose = T
  ){

  #Skip pseudocounts if not required
  if(dimsum_meta[["fitnessDropoutPseudocount"]]==0){
    return(input_dt)
  }

  if(verbose){dimsum__status_message("Adding pseudocount to output samples with dropout...\n")}

  #Add psuedocounts for all output replicates
  for(j in all_reps){
    input_dt[get(paste0("count_e", j, "_s0"))>0 & get(paste0("count_e", j, "_s1"))==0, paste0("count_e", j, "_s1") := .SD[[1]] + dimsum_meta[["fitnessDropoutPseudocount"]],,.SDcols = paste0("count_e", j, "_s1")]
  }  

  if(verbose){dimsum__status_message("Done\n")}

  return(input_dt)

}
