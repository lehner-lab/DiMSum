
#' dimsum__identify_single_aa_mutations
#'
#' Identify and annotate single AA substitutions.
#'
#' @param input_dt input data.table (required)
#' @param wt_AAseq WT amino acid sequence (required)
#'
#' @return data.table with single amino acid variants (and all silent variants)
#' @export
#' @import data.table
dimsum__identify_single_aa_mutations <- function(
  input_dt,
  wt_AAseq
  ){

  #WT AA sequences
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]

  ### Identify position and identity of single AA mutations (and all silent mutants)
  ###########################

  #Single AA mutants and all silent mutants
  singles_silent <- input_dt[Nham_aa==1 | (is.na(WT) & Nham_aa==0),]
  #Add position, mutant AA, WT AA and mean input count
  singles_silent[,Pos := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split),aa_seq]
  singles_silent[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]
  singles_silent[,WT_AA := wt_AAseq_split[Pos],aa_seq]

  #Remove unnecessary columns and rename fitness and sigma columns
  singles_silent <- singles_silent[,cbind(Pos,WT_AA,Mut,merge_seq,aa_seq,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,mean_count,.SD),,.SDcols = grep("_uncorr", names(singles_silent))]
  names(singles_silent) <- gsub("_uncorr", "", names(singles_silent))

  return(singles_silent)

}
