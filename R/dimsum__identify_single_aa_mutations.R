
#' dimsum__identify_single_aa_mutations
#'
#' Identify and annotate single AA substitutions.
#'
#' @param input_dt input data.table (required)
#' @param wt_AAseq WT amino acid sequence (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__identify_single_aa_mutations <- function(
  input_dt,
  wt_AAseq,
  report = TRUE,
  report_outpath = NULL
  ){

  #WT AA sequences
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]

  ### Identify position and identity of single AA mutations (and all silent mutants)
  ###########################

  #Single AA mutants and all silent mutants
  singles_silent <- input_dt[Nmut_aa==1 | (is.na(WT) & Nmut_aa==0),]
  #Add position, mutant AA, WT AA and mean input count
  singles_silent[,Pos := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split),aa_seq]
  singles_silent[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]
  singles_silent[,WT_AA := wt_AAseq_split[Pos],aa_seq]

  #Mean input count 
  if(report){
    d <- ggplot2::ggplot(singles_silent[Nmut_aa==1],ggplot2::aes(mean_count)) + 
      ggplot2::geom_density() + 
      ggplot2::scale_x_log10() +
      ggplot2::geom_vline(xintercept = 400)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_4_singles_meaninputcount.png"), d, width = 7, height = 5)
    #lower peak is from singles with AA changes 2 nt away
  }

  #Remove unnecessary columns and rename fitness and sigma columns
  singles_silent <- singles_silent[,cbind(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,.SD),,.SDcols = grep("_uncorr", names(singles_silent))]
  names(singles_silent) <- gsub("_uncorr", "", names(singles_silent))

  return(singles_silent)

}
