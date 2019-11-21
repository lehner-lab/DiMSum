
#' dimsum__filter_nuc_variants
#'
#' Filter out low count nucleotide variants.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt output path for plots and saved objects (required)
#' @param all_reps list of replicates to retain (required)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__filter_nuc_variants <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  message("Filtering out low count nucleotide variants...")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #Sample names
  input_samples <- names(input_dt)[grep(paste0("e[", all_reps_str, "]_s0_b.*_count$"), names(input_dt))]
  output_samples <- names(input_dt)[grep(paste0("e[", all_reps_str, "]_s1_b.*_count$"), names(input_dt))]

  #### Retain variants with input/output read counts > "fitnessMinInputCountAny" in ANY biological replicates 
  #### Retain variants with input/output read counts > "fitnessMinInputCountAll" in ALL biological replicates
  output_dt <- copy(input_dt)
  output_dt <- output_dt[rowSums(output_dt[,input_samples,with=F]>=dimsum_meta[["fitnessMinInputCountAny"]]) != 0]
  output_dt <- output_dt[rowSums(output_dt[,output_samples,with=F]>=dimsum_meta[["fitnessMinOutputCountAny"]]) != 0]
  output_dt <- output_dt[rowSums(output_dt[,input_samples,with=F]<dimsum_meta[["fitnessMinInputCountAll"]]) == 0]
  output_dt <- output_dt[rowSums(output_dt[,output_samples,with=F]<dimsum_meta[["fitnessMinOutputCountAll"]]) == 0]

  #### Set input read counts < "fitnessMinInputCountAny" to NA
  for(i in input_samples){
    output_dt[get(i) < dimsum_meta[["fitnessMinInputCountAny"]], as.character(i) := NA]
  }

  #### Set output read counts < "fitnessMinOutputCountAny" to NA
  for(i in output_samples){
    output_dt[get(i) < dimsum_meta[["fitnessMinOutputCountAny"]], as.character(i) := NA]
  }

  message("Done")

  return(output_dt)

}
