
#' dimsum__remove_internal_constant_region
#'
#' Remove internal constant region from nucleotide sequences (if same length as WT and perfect match found).
#'
#' @param input_dt input data.table (required)
#' @param wt_ccntseq case-coded WT nucleotide sequence (required)
#'
#' @return A data.table with internal constant regions removed from nucleotide sequences (if same length as WT and perfect match found)
#' @export
#' @import data.table
dimsum__remove_internal_constant_region <- function(
  input_dt,
  wt_ccntseq
  ){

  #Separate nucleotide sequences into pieces with and without constant region
  input_dt[indel==F, nt_seq_cnst := ""]
  input_dt[indel==F, nt_seq_wcnst := ""]
  idx <- 0
  for(b in strsplit(wt_ccntseq, "")[[1]]){
    idx <- idx + 1
    if(b %in% c("A", "C", "G", "T")){
      input_dt[indel==F, nt_seq_wcnst := paste0(nt_seq_wcnst, substr(nt_seq, idx, idx))]
    }else{
      input_dt[indel==F, nt_seq_cnst := paste0(nt_seq_cnst, substr(nt_seq, idx, idx))]
    }
  }

  #Indicate variants with WT internal constant region
  input_dt[indel==F, constant_region := F]
  input_dt[indel==F & nt_seq_cnst==input_dt[WT==T,nt_seq_cnst], constant_region := T]

  #Update nucleotide and amino acid sequences
  input_dt[indel==F & constant_region==T, nt_seq := nt_seq_wcnst]
  suppressWarnings(input_dt[indel==F & constant_region==T, aa_seq := as.character(Biostrings::translate(Biostrings::DNAStringSet(nt_seq)))])

  #Remove unnecessary columns
  output_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% c("nt_seq_cnst", "nt_seq_wcnst")]]

  return(output_dt)

}
