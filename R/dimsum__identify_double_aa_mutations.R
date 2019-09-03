
#' dimsum__identify_double_aa_mutations
#'
#' Identify and annotate double AA substitutions.
#'
#' @param input_dt input data.table (required)
#' @param singles_dt singles data.table (required)
#' @param wt_AAseq WT amino acid sequence (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__identify_double_aa_mutations <- function(
  input_dt,
  singles_dt,
  wt_AAseq,
  report = TRUE,
  report_outpath = NULL
  ){

  #WT AA sequences
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]

  ### Identify position and identity of double AA mutations
  ###########################

  #Double AA mutants
  doubles <- input_dt[Nmut_aa==2]
  #Add position, mutant AA, WT AA and mean input count
  doubles[,Pos1 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[1],aa_seq]
  doubles[,Pos2 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[2],aa_seq]
  doubles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
  doubles[,Mut2 := strsplit(aa_seq,"")[[1]][Pos2],aa_seq]
  doubles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
  doubles[,WT_AA2 := wt_AAseq_split[Pos2],aa_seq]
  #Mean counts
  doubles <- merge(doubles,singles_dt[,.(Pos,Mut,s1_mean_count = mean_count)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
  doubles <- merge(doubles,singles_dt[,.(Pos,Mut,s2_mean_count = mean_count)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))

  #Mean input count of doubles split accoring to mean input count of singles
  if(report){
    d <- ggplot2::ggplot(doubles,ggplot2::aes(mean_count,..count..,color = s1_mean_count > 400 & s2_mean_count > 400)) +
      ggplot2::geom_density() + 
      ggplot2::scale_x_log10() +
      ggplot2::labs(color="high single input")
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_4_doubles_meaninputcount.png"), d, width = 7, height = 5)
  }

  return(doubles)

}