
#' dimsum__identify_single_aa_mutations
#'
#' Identify and annotate single AA substitutions.
#'
#' @param input_dt input data.table (required)
#' @param wt_AAseq WT amino acid sequence (required)
#'
#' @return data.table with single amino acid variants
#' @export
#' @import data.table
dimsum__identify_single_aa_mutations <- function(
  input_dt,
  wt_AAseq
  ){

  #WT AA sequences
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]

  ### Identify position and identity of single AA mutations
  ###########################

  #Single AA mutants
  singles <- input_dt[Nham_aa==1,]

  #If no single AA mutants, return empty data.table
  if(nrow(singles)==0){
    return(data.table())
  }

  #Add position, mutant AA, WT AA and mean input count
  singles[,Pos := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split),aa_seq]
  singles[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]
  singles[,WT_AA := wt_AAseq_split[Pos],aa_seq]

  #Remove unnecessary columns and rename fitness and sigma columns
  singles <- singles[,cbind(Pos,WT_AA,Mut,merge_seq,aa_seq,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,error_model,mean_count,.SD),,.SDcols = grep("_uncorr", names(singles))]
  names(singles) <- gsub("_uncorr", "", names(singles))

  return(singles)

}
