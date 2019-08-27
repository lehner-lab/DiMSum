
#' dimsum__check_nt_fitness
#'
#' Check fitness and count-based error of single and double nucleotide substitution variants.
#'
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__check_nt_fitness <- function(
  input_dt,
  all_reps,
  report = TRUE,
  report_outpath = NULL
  ){

  message("Checking fitness and count-based error of single and double nucleotide substitution variants...")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #Reformat
  input_dt[,merge_seq := nt_seq,nt_seq]
  #Dummy column (for compatibility)
  input_dt[,Nmut_codons := 0]
  input_dt <- input_dt[,.SD,merge_seq,.SDcols = c("aa_seq","Nins_nt","Ndel_nt","Nsub_nt","Nmut_nt","Nins_aa","Ndel_aa","Nsub_aa","Nmut_aa","Nmut_codons","WT","STOP",names(input_dt)[grep(names(input_dt),pattern="_count$")])]

  #Add up counts for biological output reps
  for (E in all_reps) {
    idx <- names(input_dt)[grep(names(input_dt),pattern = paste0("e",E,"_s1_b"))]
    #Data before aggregating at AA level
    input_dt[,paste0("count_e",E,"_s1") := rowSums(.SD),,.SDcols = idx]
    names(input_dt)[grep(names(input_dt),pattern = paste0("e",E,"_s0_b"))] <- paste0("count_e",E,"_s0")
  }
  input_dt <- unique(input_dt[,.SD,merge_seq,.SDcols = c("aa_seq","Nins_nt","Ndel_nt","Nsub_nt","Nmut_nt","Nins_aa","Ndel_aa","Nsub_aa","Nmut_aa","Nmut_codons","WT","STOP",names(input_dt)[grep(names(input_dt),pattern="^count")])])

  #Calculate fitness
  input_dt <- dimsum__calculate_fitness(
    input_dt = input_dt,
    all_reps = all_reps,
    sequence_type = "noncoding",
    report = F,
    verbose = F)

  #first check potential for replicate error
  input_dt[,var_fitness := rowSums((rowMeans(.SD[,1:nchar(all_reps_str)],na.rm=T) - .SD[,1:nchar(all_reps_str)])^2/(.SD[,(nchar(all_reps_str)+1):(2*nchar(all_reps_str))]^2),na.rm=T) / 
                rowSums(1/(.SD[,(nchar(all_reps_str)+1):(2*nchar(all_reps_str))]^2 ),na.rm=T),
              ,.SDcols = c(grep(names(input_dt),pattern=paste0("fitness[", all_reps_str, "]")),grep(names(input_dt),pattern=paste0("sigma[", all_reps_str, "]")))]
  input_dt[,avg_sigma := rowMeans(.SD[,(nchar(all_reps_str)+1):(2*nchar(all_reps_str))],na.rm=T),
              ,.SDcols = c(grep(names(input_dt),pattern=paste0("fitness[", all_reps_str, "]")),grep(names(input_dt),pattern=paste0("sigma[", all_reps_str, "]")))]
  input_dt[,isNA := rowSums(is.na(.SD)),,.SDcols = grep(names(input_dt),pattern=paste0("fitness[", all_reps_str, "]"))]
  input_dt[isNA==0 & var_fitness != Inf & avg_sigma != Inf,avg_sigma_fit := loess(var_fitness ~ avg_sigma,span=0.75)$fitted]

  replicate_error <- input_dt[,min(sqrt(avg_sigma_fit),na.rm=T)]
  message(paste0("Replicate error is:", replicate_error))

  #Average sigma versus fitness replicate error (hexbin for single and double nucleotide variants only)
  if(report){
    d <- ggplot2::ggplot(input_dt[isNA==0 & !is.infinite(avg_sigma)],ggplot2::aes(avg_sigma)) +
      ggplot2::stat_binhex(data = input_dt[isNA==0 & !is.infinite(avg_sigma) & between(Nmut_nt, 1, 2)], ggplot2::aes(y=sqrt(var_fitness), color = as.factor(Nmut_nt)), bins=100, size = 0.2) +
      ggplot2::labs(color = "Nmut_nt") +
      ggplot2::scale_fill_gradientn(colours=c("white", "black")) +
      # ggplot2::geom_point(ggplot2::aes(y=sqrt(var_fitness))) +
      ggplot2::geom_line(ggplot2::aes(y=sqrt(avg_sigma_fit)),color="black",linetype=2) +
      ggplot2::geom_line(ggplot2::aes(y=sqrt(avg_sigma^2 + replicate_error^2)),color="black") +
      ggplot2::geom_abline(color="darkgrey",linetype=2) + 
      ggplot2::scale_x_log10() + 
      ggplot2::scale_y_log10() +
      ggplot2::xlab("Average count-based (Poisson) fitness error") + ggplot2::ylab("Fitness standard deviation between biological replicates") +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor=ggplot2::element_blank())
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_singledouble_fitness_replicateerror_vs_avgsigma.png"), d, width = 6, height = 5)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_singledouble_fitness_replicateerror_vs_avgsigma.pdf"), d, width = 6, height = 5)

    #Add pseudocount of 10^-3 to sqrt(var_fitness)
    d <- ggplot2::ggplot(input_dt[isNA==0 & !is.infinite(avg_sigma)],ggplot2::aes(avg_sigma)) +
      ggplot2::stat_binhex(data = input_dt[isNA==0 & !is.infinite(avg_sigma) & between(Nmut_nt, 1, 2)], ggplot2::aes(y=sqrt(var_fitness)+10^-3, color = as.factor(Nmut_nt)), bins=100, size = 0.1) +
      ggplot2::labs(color = "Nmut_nt") +
      ggplot2::scale_fill_gradientn(colours=c("white", "black")) +
      # ggplot2::geom_point(ggplot2::aes(y=sqrt(var_fitness))) +
      ggplot2::geom_line(ggplot2::aes(y=sqrt(avg_sigma_fit)),color="black",linetype=2) +
      ggplot2::geom_line(ggplot2::aes(y=sqrt(avg_sigma^2 + replicate_error^2)),color="black") +
      ggplot2::geom_abline(color="darkgrey",linetype=2) + 
      ggplot2::scale_x_log10() + 
      ggplot2::scale_y_log10() +
      ggplot2::xlab("Average count-based (Poisson) fitness error") + ggplot2::ylab("Fitness standard deviation between biological replicates") +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor=ggplot2::element_blank())
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_singledouble_fitness_replicateerror_vs_avgsigma_pseudocount.png"), d, width = 6, height = 5)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_1_singledouble_fitness_replicateerror_vs_avgsigma_pseudocount.pdf"), d, width = 6, height = 5)
  }

  message("Done")

}
