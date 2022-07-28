
#' dimsum__identify_single_nt_mutations
#'
#' Identify and annotate single nucleotide substitutions.
#'
#' @param input_dt input data.table (required)
#' @param wt_ntseq WT nucleotide sequence (required)
#'
#' @return data.table with single nucleotide variants
#' @export
#' @import data.table
dimsum__identify_single_nt_mutations <- function(
  input_dt,
  wt_ntseq
  ){

  #WT nucleotide sequences
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]

  ### Identify position and identity of single nucleotide mutations
  ###########################

  #Single nucleotide mutants
  singles_silent <- input_dt[Nham_nt==1,]

  #If no single nucleotide mutants, return empty data.table
  if(nrow(singles_silent)==0){
    return(data.table())
  }

  #Add position, mutant nucleotide, WT nucleotide (WT_AA)
  singles_silent[,Pos := which(strsplit(merge_seq,"")[[1]] !=wt_ntseq_split),merge_seq]
  singles_silent[,Mut := strsplit(merge_seq,"")[[1]][Pos],merge_seq]
  singles_silent[,WT_AA := wt_ntseq_split[Pos],merge_seq]

  #Remove unnecessary columns and rename fitness and sigma columns
  singles_silent <- singles_silent[,cbind(Pos,WT_AA,Mut,merge_seq,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,error_model,mean_count,.SD),,.SDcols = grep("_uncorr", names(singles_silent))]
  names(singles_silent) <- gsub("_uncorr", "", names(singles_silent))

  return(singles_silent)

}
