
#' dimsum__aggregate_AA_variants
#'
#' Aggregate counts from variants that are identical at the AA level and without synonymous mutations.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#'
#' @return data.table with aggregated counts
#' @export
#' @import data.table
dimsum__aggregate_AA_variants <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  dimsum__status_message("Aggregating counts for identical amino acid variants...\n")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #Sample names
  input_samples <- names(input_dt)[grep(paste0("e[", all_reps_str, "]_s0_b.*_count$"), names(input_dt))]
  output_samples <- names(input_dt)[grep(paste0("e[", all_reps_str, "]_s1_b.*_count$"), names(input_dt))]

  dimsum__check_variants(dimsum_meta = dimsum_meta, input_dt = input_dt)

  #Set merge_seq to nt_seq
  input_dt[,merge_seq := aa_seq,nt_seq]

  #Define silent/synonymous variants as WT
  input_dt[Nham_aa==0,WT := T]
  input_dt[Nham_aa==0,Nham_nt := 0]

  #For all count columns
  for(j in c(input_samples, output_samples)){
    #Aggregate counts accross identical AA variants
    input_dt[!is.na(get(j)),paste0(j,"_agg") := sum(.SD, na.rm = T),merge_seq,.SDcols = j]
    #Set all synonyms to the same count
    input_dt[,paste0(j,"_agg") := unique(na.omit(.SD)),merge_seq,.SDcols = paste0(j,"_agg")]
  }

  #Retain only one row per AA variant
  output_dt <- input_dt[!duplicated(merge_seq),.SD,merge_seq,.SDcols = c(
    "nt_seq","aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","indel","STOP","STOP_readthrough",
    names(input_dt)[grep(names(input_dt),pattern="_agg$")])]

  #Revert to original names of aggregated count columns
  names(output_dt)[grep(names(output_dt),pattern="_agg$")] <- gsub("_agg$", "", names(output_dt)[grep(names(output_dt),pattern="_agg$")])

  #Nham_nt and Nmut_codons columns are meaningless for nonsynonymous variants after aggregation
  output_dt[Nham_aa!=0, Nham_nt := NA]
  output_dt[Nham_aa!=0, Nmut_codons := NA]

  #nt_seq column is meaningless after aggregation
  output_dt[, nt_seq := NA]

  dimsum__status_message("Done\n")

  return(output_dt)

}
