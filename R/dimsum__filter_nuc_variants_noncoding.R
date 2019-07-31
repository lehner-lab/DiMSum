
#' dimsum__filter_nuc_variants_noncoding
#'
#' Filter nucleotide variants from noncoding sequences.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt output path for plots and saved objects (required)
#' @param wt_ntseq WT nucleotide sequence (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__filter_nuc_variants_noncoding <- function(
  dimsum_meta,
  input_dt,
  wt_ntseq,
  report = TRUE,
  report_outpath = NULL
  ){

  message("Filtering low count nucleotide variants...")

  #WT nucleotide sequences
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]

  #Sample names
  input_samples <- names(input_dt)[grep("_e.*_s0_b.*_count$", names(input_dt))]

  ### Minimum input read count in any replicate (to remove very low confidence variants)
  ###########################

  #Minimum input read count
  min_input_read_count <- 5

  ### Retain nucleotide variants with same length as WT and max 2 nucleotide substitutions only
  ### Retain only variants with >=min_input_read_count input readcounts in any biological replicate 
  ###########################

  #Retain nucleotide variants with same length as WT and max 2 nucleotide substitutions only
  input_dt <- input_dt[nchar(nt_seq)==nchar(input_dt[WT==T,nt_seq]) & Nmut_nt<=2 & Nmut_nt==Nsub_nt,]

  #Add number of codons affected by mutations
  input_dt[,Nmut_codons := length(unique(ceiling(which(strsplit(nt_seq,"")[[1]] != wt_ntseq_split)/3))),nt_seq]

  #Example histogram of input1 counts split by number of nucleotide mutations 
  if(report){
    temp_input_name <- names(input_dt)[grep("_e1_s0_b.*_count$", names(input_dt))][1]
    d <- ggplot2::ggplot(input_dt[between(Nmut_nt,1,4) & get(temp_input_name) > 0],
           ggplot2::aes(get(temp_input_name),color=factor(Nmut_nt))) +
      ggplot2::geom_density() +
      ggplot2::scale_x_log10()
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_input1_count_hist.png"), d, width = 7, height = 5)
  }

  #### only retain variants with at least 5 input readcounts in any of the biological replicates
  output_dt <- copy(input_dt[rowSums(input_dt[,input_samples,with=F]>min_input_read_count) != 0])

  #Plot pairwise input sample count correlations for all single mutants (non-synonymous only)
  if(report){
    set.seed(1)
    d <- GGally::ggpairs(log10(output_dt[Nmut_nt==1][sample(x = .N,min(c(.N,10000))),grep(names(output_dt),pattern="input"),with=F]+1),
            upper=list(continuous = "cor"))
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_2_scatterplotmatrix_singlemutants.png"), d, width = 10, height = 10)
  }
  
  #Plot pairwise input sample count correlations for random sample of 10k variants
  if(report){
    temp <- output_dt[is.na(WT),][sample(x = .N,min(c(.N,10000)))]
    d <- GGally::ggpairs(cbind(log10(temp[,grep(names(output_dt),pattern="input"),with=F]+1), Nmut_nt=as.factor(temp[,Nmut_nt])),
      columns = 1:length(grep(names(output_dt),pattern="input")),
      mapping = ggplot2::aes(color = Nmut_nt),
      upper=list(continuous = "cor"))
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_2_scatterplotmatrix_random10kmutants.png"), d, width = 10, height = 10)
  }

  message("Done")

  return(output_dt)

}
