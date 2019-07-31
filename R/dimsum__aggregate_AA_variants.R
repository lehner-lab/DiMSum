
#' dimsum__aggregate_AA_variants
#'
#' Aggregate counts from variants that are identical at the AA level and without synonymous mutations.
#'
#' @param input_dt output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__aggregate_AA_variants <- function(
  input_dt
  ){

  message("Aggregating counts for identical amino acid variants...")

  #Add merge_seq (aa_seq if nonsynonymous mutations, otherwise nt_seq)
  input_dt[,merge_seq := aa_seq,nt_seq]
  input_dt[Nmut_aa==0,merge_seq := nt_seq,nt_seq]

  #For all count columns
  idx <- names(input_dt)[grep(names(input_dt),pattern="_count$")]
  for (i in seq_along(idx)) {
    #Aggregate counts accross identical AA variants
    input_dt[,paste0(idx[i],"_agg") := sum(.SD),merge_seq,.SDcols = idx[i]]
  }
  #Retain only one row per AA variant
  output_dt <- input_dt[!duplicated(merge_seq),.SD,merge_seq,.SDcols = c("aa_seq","Nmut_nt","Nmut_aa","Nmut_codons","WT","STOP",names(input_dt)[grep(names(input_dt),pattern="_agg$")])]

  message("Done")

  return(output_dt)

}
