
#' dimsum__filter_nuc_variants
#'
#' Filter out low count nucleotide variants.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#'
#' @return data.table with low count nucleotide variants filtered out
#' @export
#' @import data.table
dimsum__filter_nuc_variants <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  dimsum__status_message("Filtering out low count nucleotide variants...\n")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Sample names
  input_samples <- names(input_dt)[grep(paste0("e(", all_reps_str, ")_s0_b.*_count$"), names(input_dt))]
  output_samples <- names(input_dt)[grep(paste0("e(", all_reps_str, ")_s1_b.*_count$"), names(input_dt))]

  #### Retain variants with input/output read counts > "fitnessMinInputCountAny" in ANY biological replicates 
  #### Retain variants with input/output read counts > "fitnessMinInputCountAll" in ALL biological replicates
  output_dt <- copy(input_dt)
  #fitnessMinInputCountAny
  if(is.list(dimsum_meta[["fitnessMinInputCountAny"]])){
    for(i in names(dimsum_meta[["fitnessMinInputCountAny"]])){
      output_dt <- output_dt[(Nham_nt == as.integer(i) & rowSums(output_dt[,input_samples,with=F]>=dimsum_meta[["fitnessMinInputCountAny"]][[i]]) != 0) | Nham_nt != as.integer(i)]
    }
  }else{
    output_dt <- output_dt[rowSums(output_dt[,input_samples,with=F]>=dimsum_meta[["fitnessMinInputCountAny"]]) != 0]
  }
  #fitnessMinOutputCountAny
  if(is.list(dimsum_meta[["fitnessMinOutputCountAny"]])){
    for(i in names(dimsum_meta[["fitnessMinOutputCountAny"]])){
      output_dt <- output_dt[(Nham_nt == as.integer(i) & rowSums(output_dt[,output_samples,with=F]>=dimsum_meta[["fitnessMinOutputCountAny"]][[i]]) != 0) | Nham_nt != as.integer(i)]
    }
  }else{
    output_dt <- output_dt[rowSums(output_dt[,output_samples,with=F]>=dimsum_meta[["fitnessMinOutputCountAny"]]) != 0]
  }
  #fitnessMinInputCountAll
  if(is.list(dimsum_meta[["fitnessMinInputCountAll"]])){
    for(i in names(dimsum_meta[["fitnessMinInputCountAll"]])){
      output_dt <- output_dt[(Nham_nt == as.integer(i) & rowSums(output_dt[,input_samples,with=F]<dimsum_meta[["fitnessMinInputCountAll"]][[i]]) == 0) | Nham_nt != as.integer(i)]
    }
  }else{
    output_dt <- output_dt[rowSums(output_dt[,input_samples,with=F]<dimsum_meta[["fitnessMinInputCountAll"]]) == 0]
  }
  #fitnessMinOutputCountAll
  if(is.list(dimsum_meta[["fitnessMinOutputCountAll"]])){
    for(i in names(dimsum_meta[["fitnessMinOutputCountAll"]])){
      output_dt <- output_dt[(Nham_nt == as.integer(i) & rowSums(output_dt[,output_samples,with=F]<dimsum_meta[["fitnessMinOutputCountAll"]][[i]]) == 0) | Nham_nt != as.integer(i)]
    }
  }else{
    output_dt <- output_dt[rowSums(output_dt[,output_samples,with=F]<dimsum_meta[["fitnessMinOutputCountAll"]]) == 0]
  }

  #### Set input read counts < "fitnessMinInputCountAny" to NA
  for(i in input_samples){
    if(is.list(dimsum_meta[["fitnessMinInputCountAny"]])){
      for(j in names(dimsum_meta[["fitnessMinInputCountAny"]])){
        output_dt[Nham_nt == as.integer(j) & get(i) < dimsum_meta[["fitnessMinInputCountAny"]][[j]], as.character(i) := NA]
      }
    }else{
      output_dt[get(i) < dimsum_meta[["fitnessMinInputCountAny"]], as.character(i) := NA]
    }
  }

  #### Set output read counts < "fitnessMinOutputCountAny" to NA
  for(i in output_samples){
    if(is.list(dimsum_meta[["fitnessMinOutputCountAny"]])){
      for(j in names(dimsum_meta[["fitnessMinOutputCountAny"]])){
        output_dt[Nham_nt == as.integer(j) & get(i) < dimsum_meta[["fitnessMinOutputCountAny"]][[j]], as.character(i) := NA]
      }
    }else{
      output_dt[get(i) < dimsum_meta[["fitnessMinOutputCountAny"]], as.character(i) := NA]
    }
  }

  dimsum__check_variants(dimsum_meta = dimsum_meta, input_dt = output_dt)

  dimsum__status_message("Done\n")

  return(output_dt)

}
