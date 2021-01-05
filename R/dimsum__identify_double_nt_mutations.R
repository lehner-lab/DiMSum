
#' dimsum__identify_double_nt_mutations
#'
#' Identify and annotate double nucleotide substitutions.
#'
#' @param input_dt input data.table (required)
#' @param singles_dt singles data.table (required)
#' @param wt_ntseq WT nucleotide sequence (required)
#'
#' @return data.table with double nucleotide variants
#' @export
#' @import data.table
dimsum__identify_double_nt_mutations <- function(
  input_dt,
  singles_dt,
  wt_ntseq
  ){

  #WT nucleotide sequences
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]

  ### Identify position and identity of double nucleotide mutations
  ###########################

  #Double nucleotide mutants
  doubles <- input_dt[Nham_nt==2]

  #If no double nucleotide mutants, return empty data.table
  if(nrow(doubles)==0){
    return(data.table())
  }

  #Add positions, mutant nucleotides, WT nucleotides (WT_AA1, WT_AA2) and mean input count
  doubles[,Pos1 := which(strsplit(merge_seq,"")[[1]] !=wt_ntseq_split)[1],merge_seq]
  doubles[,Pos2 := which(strsplit(merge_seq,"")[[1]] !=wt_ntseq_split)[2],merge_seq]
  doubles[,Mut1 := strsplit(merge_seq,"")[[1]][Pos1],merge_seq]
  doubles[,Mut2 := strsplit(merge_seq,"")[[1]][Pos2],merge_seq]
  doubles[,WT_AA1 := wt_ntseq_split[Pos1],merge_seq]
  doubles[,WT_AA2 := wt_ntseq_split[Pos2],merge_seq]

  #Mean counts
  if(nrow(singles_dt)!=0){
    doubles <- merge(doubles,singles_dt[,.(Pos,Mut,s1_mean_count = mean_count)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    doubles <- merge(doubles,singles_dt[,.(Pos,Mut,s2_mean_count = mean_count)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
  }else{
    doubles[, s1_mean_count := NA]
    doubles[, s2_mean_count := NA]
  }

  return(doubles)

}
