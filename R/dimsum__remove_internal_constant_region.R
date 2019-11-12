
#' dimsum__remove_internal_constant_region
#'
#' Remove internal constant region from nucleotide sequences (if perfect match found and otherwise removed entirely).
#'
#' @param input_dt input data.table (required)
#' @param wt_ccntseq case-coded WT nucleotide sequence (required)
#'
#' @return A data.table with internal constant regions removed from nucleotide sequences (if perfect match found and otherwise removed entirely)
#' @export
#' @import data.table
dimsum__remove_internal_constant_region <- function(
  input_dt,
  wt_ccntseq
  ){

  message("Removing internal constant region from nucleotide variants...")

  #Separate nucleotide sequences into pieces with and without constant region
  input_dt[, nt_seq_cnst := ""]
  input_dt[, nt_seq_wcnst := ""]
  idx <- 0
  for(b in strsplit(wt_ccntseq, "")[[1]]){
    idx <- idx + 1
    if(b %in% c("A", "C", "G", "T")){
      input_dt[, nt_seq_wcnst := paste0(nt_seq_wcnst, substr(nt_seq, idx, idx))]
    }else{
      input_dt[, nt_seq_cnst := paste0(nt_seq_cnst, substr(nt_seq, idx, idx))]
    }
  }

  #Subset to variants with WT internal constant region
  input_dt <- input_dt[nt_seq_cnst==input_dt[WT==T,nt_seq_cnst],]

  #Update nucleotide and amino acid sequences
  input_dt[, nt_seq := nt_seq_wcnst]
  suppressWarnings(input_dt[, aa_seq := as.character(Biostrings::translate(Biostrings::DNAStringSet(nt_seq)))])

  #Indicate STOPs
  input_dt[,STOP := ifelse(length(grep(aa_seq,pattern="\\*"))==1,TRUE,FALSE),aa_seq]

  #Remove unnecessary columns
  output_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% c("nt_seq_cnst", "nt_seq_wcnst")]]

  message("Done")

  return(output_dt)

}
