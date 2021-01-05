
#' dimsum__identify_double_aa_mutations
#'
#' Identify and annotate double AA substitutions.
#'
#' @param input_dt input data.table (required)
#' @param singles_dt singles data.table (required)
#' @param wt_AAseq WT amino acid sequence (required)
#'
#' @return data.table with double amino acid variants
#' @export
#' @import data.table
dimsum__identify_double_aa_mutations <- function(
  input_dt,
  singles_dt,
  wt_AAseq
  ){

  #WT AA sequences
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]

  ### Identify position and identity of double AA mutations
  ###########################

  #Double AA mutants
  doubles <- input_dt[Nham_aa==2]

  #If no double AA mutants, return empty data.table
  if(nrow(doubles)==0){
    return(data.table())
  }

  #Add position, mutant AA, WT AA and mean input count
  doubles[,Pos1 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[1],aa_seq]
  doubles[,Pos2 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[2],aa_seq]
  doubles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
  doubles[,Mut2 := strsplit(aa_seq,"")[[1]][Pos2],aa_seq]
  doubles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
  doubles[,WT_AA2 := wt_AAseq_split[Pos2],aa_seq]

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
