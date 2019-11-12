
#' dimsum__filter_variants_coding
#'
#' Filter for desired coding variants.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt output path for plots and saved objects (required)
#' @param wt_ntseq WT nucleotide sequence (required)
#' @param all_reps list of replicates to retain (required)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__filter_variants_coding <- function(
  dimsum_meta,
  input_dt,
  wt_ntseq,
  all_reps){

  message("Filtering for desired coding variants...")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #WT nucleotide sequences
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]

  #Sample names
  input_samples <- names(input_dt)[grep(paste0("e[", all_reps_str, "]_s0_b.*_count$"), names(input_dt))]

  ### Retain nucleotide variants with max dimsum_meta[["fitnessMaxSubstitutions"]] amino acid mutations only
  ###########################

  #Retain nucleotide variants with max dimsum_meta[["fitnessMaxSubstitutions"]] amino acid mutations only
  input_dt <- input_dt[Nham_aa<=dimsum_meta[["fitnessMaxSubstitutions"]],]

  #Add number of codons affected by mutations
  input_dt[,Nmut_codons := length(unique(ceiling(which(strsplit(nt_seq,"")[[1]] != wt_ntseq_split)/3))),nt_seq]

  ### Output data.table
  ### Retain only either purely nonsynonymous mutations or purely silent mutations
  output_dt <- copy(input_dt[((Nmut_codons-Nham_aa) == 0 | Nham_aa == 0)])
  
  message("Done")

  return(output_dt)

}
