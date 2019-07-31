
#' dimsum__detect_sequence_type
#'
#' Identify and annotate double AA substitutions.
#'
#' @param input_sequence input nucleotide sequence (required)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__detect_sequence_type <- function(
  input_sequence
  ){

  #Sequence length a multiple of 3?
  if(nchar(input_sequence) %% 3 == 0){
    aa_seq_split <- seqinr::translate(strsplit(input_sequence,split="")[[1]])
    #Translated sequence contains STOP?
    if(!"*" %in% aa_seq_split){
      return("coding")
    }
  }
  
  #Sequence length not a multiple of 3 or contains STOP
  return("noncoding")

}
