
#' dimsum__aggregate_AA_variants_fitness
#'
#' Aggregate fitness and error from variants that are identical at the AA level and without synonymous mutations.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#'
#' @return data.table with aggregated counts
#' @export
#' @import data.table
dimsum__aggregate_AA_variants_fitness <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  message("Aggregating fitness for identical amino acid variants...")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #Sample names
  input_samples <- names(input_dt)[grep(paste0("count_e[", all_reps_str, "]_s0$"), names(input_dt))]
  output_samples <- names(input_dt)[grep(paste0("count_e[", all_reps_str, "]_s1$"), names(input_dt))]

  dimsum__check_variants(dimsum_meta = dimsum_meta, input_dt = input_dt)

  #Set merge_seq to nt_seq
  input_dt[,merge_seq := aa_seq,nt_seq]

  #Set merge_seq to nt_seq if WT or silent/synonymous variants
  input_dt[Nham_aa==0,merge_seq := nt_seq,nt_seq]

  #For all count columns
  for(j in c(input_samples, output_samples)){
    #Aggregate counts accross identical AA variants
    input_dt[!is.na(get(j)),paste0(j,"_agg") := sum(.SD, na.rm = T),merge_seq,.SDcols = j]
  }

  #Calculate mean input counts
  input_dt[,mean_count := rowMeans(.SD, na.rm = T),,.SDcols = paste0("count_e", all_reps, "_s0_agg")]

  #For all fitness columns
  idx <- names(input_dt)[grep(names(input_dt),pattern="^fitness")]
  for (i in seq_along(idx)) {
    #Aggregate fitness accross identical AA variants
    input_dt[!is.na(get(idx[i])) & !is.na(get(gsub("fitness", "sigma", idx[i]))),paste0(idx[i],"_agg") := sum(.SD[[1]]/(.SD[[2]]^2), na.rm = T)/sum(1/(.SD[[2]]^2), na.rm = T),merge_seq,.SDcols = c(idx[i], gsub("fitness", "sigma", idx[i]))]
    input_dt[!is.na(get(gsub("fitness", "sigma", idx[i]))),paste0(gsub("fitness", "sigma", idx[i]),"_agg") := sqrt(1/sum(1/(.SD[[1]]^2), na.rm = T)),merge_seq,.SDcols = gsub("fitness", "sigma", idx[i])]
  }

  #Retain only one row per AA variant
  output_dt <- input_dt[!duplicated(merge_seq),.SD,merge_seq,.SDcols = c(
    "aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","STOP","STOP_readthrough",
    "mean_count",
    names(input_dt)[grep(names(input_dt),pattern="_agg$")])]

  #Revert to original names of aggregated count columns
  names(output_dt)[grep(names(output_dt),pattern="_agg$")] <- gsub("_agg$", "", names(output_dt)[grep(names(output_dt),pattern="_agg$")])

  message("Done")

  return(output_dt)

}
