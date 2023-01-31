
#' dimsum__filter_variants_noncoding
#'
#' Filter for desired non-coding variants.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param wt_ntseq WT nucleotide sequence (required)
#' @param all_reps list of replicates to retain (required)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__filter_variants_noncoding <- function(
  dimsum_meta,
  input_dt,
  wt_ntseq,
  all_reps){

  # dimsum__status_message("Filtering for desired non-coding variants...\n")

  # #Number of input and output replicates
  # all_reps_str <- paste0(all_reps, collapse="|")

  # #WT nucleotide sequences
  # wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]

  # #Sample names
  # input_samples <- names(input_dt)[grep(paste0("e(", all_reps_str, ")_s0_b.*_count$"), names(input_dt))]

  ### Retain nucleotide variants with max dimsum_meta[["maxSubstitutions"]] nucleotide mutations only
  ###########################

  #Retain nucleotide variants with max dimsum_meta[["maxSubstitutions"]] nucleotide mutations only
  input_dt <- input_dt[Nham_nt<=dimsum_meta[["maxSubstitutions"]],]

  # #Add number of codons affected by mutations
  # input_dt[,Nmut_codons := length(unique(ceiling(which(strsplit(nt_seq,"")[[1]] != wt_ntseq_split)/3))),nt_seq]

  # ### Output data.table
  # output_dt <- copy(input_dt)
  
  # dimsum__status_message("Done\n")

  # return(output_dt)

}
