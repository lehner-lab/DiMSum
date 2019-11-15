
#' dimsum__identify_STOP_mutations
#'
#' Identify and annotate single nucleotide substitutions.
#'
#' @param input_dt input data.table (required)
#'
#' @return data.table with single nucleotide variants
#' @export
#' @import data.table
dimsum__identify_STOP_mutations <- function(
  input_dt
  ){

  #WT AA sequence
  wt_AAseq <- input_dt[WT==T,aa_seq]

  #If STOPs in WT remove first
  if(grepl("\\*", wt_AAseq)){
    input_dt[, aa_seq_ns := ""]
    idx <- 0
    for(b in strsplit(wt_AAseq, "")[[1]]){
      idx <- idx + 1
      if(b!="*"){input_dt[, aa_seq_ns := paste0(aa_seq_ns, substr(aa_seq, idx, idx))]}
    }
  }else{
    input_dt[, aa_seq_ns := aa_seq]
  }
  #Indicate STOPs not present in WT
  input_dt[,STOP := ifelse(length(grep(aa_seq_ns,pattern="\\*"))==1,TRUE,FALSE),aa_seq_ns]

  #Indicate readthrough STOP mutations (if terminating STOP exists in WT)
  input_dt[, STOP_readthrough := FALSE]
  if(input_dt[WT==T,substr(aa_seq, nchar(wt_AAseq), nchar(wt_AAseq))=="*"]){
    input_dt[, STOP_readthrough := (substr(aa_seq, nchar(wt_AAseq), nchar(wt_AAseq))!="*")]
  }

  return(input_dt[,.SD,,.SDcols = names(input_dt)[names(input_dt)!="aa_seq_ns"]])

}
