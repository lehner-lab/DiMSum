
#' dimsum__identify_STOP_mutations
#'
#' Identify and annotate STOP codon variants.
#'
#' @param input_dt input data.table (required)
#'
#' @return data.table with STOP codon variants
#' @export
#' @import data.table
dimsum__identify_STOP_mutations <- function(
  input_dt
  ){

  #WT AA sequence
  wt_AAseq <- input_dt[WT==T,aa_seq]

  #If WT has terminating STOP remove last residue in all sequences
  if(substr(wt_AAseq, nchar(wt_AAseq), nchar(wt_AAseq))=="*"){
    input_dt[, aa_seq_ns := substr(aa_seq, 1, nchar(aa_seq)-1)]
  }else{
    input_dt[, aa_seq_ns := aa_seq]
  }
  #Indicate premature STOP (regardless of presence in WT)
  input_dt[,STOP := ifelse(length(grep(aa_seq_ns,pattern="\\*"))==1,TRUE,FALSE),aa_seq_ns]

  #Indicate readthrough STOP mutations (if terminating STOP exists in WT)
  input_dt[, STOP_readthrough := FALSE]
  if(substr(wt_AAseq, nchar(wt_AAseq), nchar(wt_AAseq))=="*"){
    input_dt[, STOP_readthrough := (substr(aa_seq, nchar(aa_seq), nchar(aa_seq))!="*")]
  }

  return(input_dt[,.SD,,.SDcols = names(input_dt)[names(input_dt)!="aa_seq_ns"]])

}
