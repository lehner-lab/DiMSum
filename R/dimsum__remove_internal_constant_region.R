
#' dimsum__remove_internal_constant_region
#'
#' Remove internal constant region from nucleotide sequences (if perfect match found).
#'
#' @param input_dt output path for plots and saved objects (required)
#' @param wt_ccntseq case-coded WT nucleotide sequence (required)
#'
#' @return A data.table with internal constant regions removed and 
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

  #Remove unnecessary columns
  output_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% c("nt_seq_cnst", "nt_seq_wcnst")]]

  #WT nucleotide sequences
  wt_ntseq <- output_dt[WT==T,nt_seq]

  #WT AA sequences
  wt_AAseq <- output_dt[WT==T,aa_seq]

  #Calculate number of aa mutations (insertions, deletions, substitutions)
  mut_counts <- attr(utils::adist(output_dt[,aa_seq], wt_AAseq, counts = T), "counts")
  output_dt[,Nins_aa := mut_counts[,1,2]]
  output_dt[,Ndel_aa := mut_counts[,1,1]]
  output_dt[,Nsub_aa := mut_counts[,1,3]]
  output_dt[,Nmut_aa := Nins_aa+Ndel_aa+Nsub_aa]
  #Calculate number of nucleotide mutations (insertions, deletions, substitutions)
  mut_counts <- attr(utils::adist(output_dt[,nt_seq], wt_ntseq, counts = T), "counts")
  output_dt[,Nins_nt := mut_counts[,1,2]]
  output_dt[,Ndel_nt := mut_counts[,1,1]]
  output_dt[,Nsub_nt := mut_counts[,1,3]]
  output_dt[,Nmut_nt := Nins_nt+Ndel_nt+Nsub_nt]

  message("Done")

  return(output_dt)

}
