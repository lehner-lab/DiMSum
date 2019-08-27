
#' dimsum__calculate_fitness
#'
#' Calculate fitness and count-based error.
#'
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param sequence_type coding potential of sequence: noncoding/coding (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#' @param verbose whether or not to print status messages (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__calculate_fitness <- function(
  input_dt,
  all_reps,
  sequence_type,
  report = TRUE,
  report_outpath = NULL,
  verbose = T
  ){

  if(verbose){message("Calculating fitness and count-based error...")}

  #For all input replicates
  for (E in all_reps) {
    f_wt_corr <- input_dt[WT == T,
                            log(.SD[,1]/.SD[,2]),,
                            .SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    input_dt[,paste0("fitness",E,"_uncorr") := log(unlist(.SD[,1]) / unlist(.SD[,2])) - rep(t(f_wt_corr),nrow(input_dt)),
            ,.SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    s_wt_counts <- input_dt[WT == T,
                              1/.SD[,1] + 1/.SD[,2],,
                              .SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    input_dt[,paste0("sigma",E,"_uncorr") := sqrt(1/unlist(.SD[,1]) + 1/unlist(.SD[,2]) + rep(t(s_wt_counts),nrow(input_dt))),
            ,.SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
  }
  #Remove unnecessary columns
  output_dt <- input_dt[,.SD,merge_seq,.SDcols = c("aa_seq","Nins_nt","Ndel_nt","Nsub_nt","Nmut_nt","Nins_aa","Ndel_aa","Nsub_aa","Nmut_aa","Nmut_codons","WT","STOP",names(input_dt)[grep(names(input_dt),pattern="^count|^fitness|^sigma")])]

  if(report){

    #Result table with fitness and count-based error for all input replicates
    DT <- NULL
    for (E in all_reps) {
      DT <- rbind(DT, output_dt[,.(count_in = get(paste0("count_e", E, "_s0")),sigma = get(paste0("sigma", E, "_uncorr")),fitness = get(paste0("fitness", E, "_uncorr")),rep = paste0("rep", E),Nmut_aa,Nmut_nt)])
    }

    #Fitness density plots
    if(sequence_type=="coding"){
      d <- ggplot2::ggplot(DT[Nmut_aa != 0 & !is.infinite(fitness) & !is.na(fitness)],ggplot2::aes(fitness,color=factor(Nmut_aa),linetype = count_in > 10))
    }else{
      d <- ggplot2::ggplot(DT[Nmut_nt != 0 & !is.infinite(fitness) & !is.na(fitness)],ggplot2::aes(fitness,color=factor(Nmut_nt),linetype = count_in > 10))
    }
    d <- d + ggplot2::geom_density()
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_3_fitness_densities.png"), d, width = 7, height = 5)

    #Fitness vs. input count hexagonal heatmap
    d <- ggplot2::ggplot(DT[!is.infinite(fitness) & !is.na(fitness)],ggplot2::aes(count_in,fitness)) +
      ggplot2::geom_hex() +
      ggplot2::scale_x_log10() +
      ggplot2::facet_wrap( ~ rep)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_3_fitness_count_hexbin.png"), d, width = 7, height = 5)

    #Sigma vs. fitness  hexagonal heatmap
    d <- ggplot2::ggplot(DT[!is.infinite(fitness) & !is.na(fitness)],ggplot2::aes(fitness,sigma)) +
      ggplot2::geom_hex() +
      ggplot2::scale_y_log10() +
      ggplot2::facet_wrap( ~ rep)
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_3_sigma_fitness_hexbin.png"), d, width = 7, height = 5)
  
  }

  #Calculate mean input counts
  output_dt[,mean_count := rowMeans(.SD),,.SDcols = paste0("count_e", all_reps, "_s0")]

  if(verbose){message("Done")}

  return(output_dt)

}
