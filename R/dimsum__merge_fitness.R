
#' dimsum__merge_fitness
#'
#' Calculate fitness and count-based error.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param doubles_dt doubles data.table (required)
#' @param singles_dt singles data.table (required)
#' @param wt_dt WT data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param min_mean_input_read_count minimum mean input read count for high confidence variants (required)
#' @param bayesian_double_fitness whether Bayesian double mutant fitness estimates exist (required)
#' @param fitness_outpath output path for saved objects (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__merge_fitness <- function(
  dimsum_meta,
  input_dt,
  doubles_dt,
  singles_dt,
  wt_dt,
  all_reps,
  min_mean_input_read_count,
  bayesian_double_fitness,
  fitness_outpath,
  report = TRUE,
  report_outpath = NULL
  ){

  message("Merging fitness estimates from biological replicates...")

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

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
      ggplot2::stat_binhex(data = input_dt[isNA==0 & !is.infinite(avg_sigma) & between(Nmut_nt, 1, 2)], ggplot2::aes(y=sqrt(var_fitness), color = as.factor(Nmut_nt)), bins=100, size = 0.1) +
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
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_5_fitness_replicateerror_vs_avgsigma.png"), d, width = 6, height = 5)
  }

  #### singles
  fitness_rx <- singles_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]"),colnames(singles_dt))]
  sigma_rx <- sqrt(singles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]"),colnames(singles_dt))]^2 + 
                    matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
  # sigma_s2 <- random_effect_model(fitness_rx,sigma_rx)
  # singles_dt[,fitness := rowSums(fitness_rx/(sigma_rx^2 + sigma_s2),na.rm=T)/rowSums(1/(sigma_rx^2 + sigma_s2),na.rm=T)]
  singles_dt[,fitness := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  # singles_dt[,sigma := sqrt(1/rowSums(1/(sigma_rx^2+sigma_s2),na.rm=T))]
  singles_dt[,sigma := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  if(report){
    d <- ggplot2::ggplot(singles_dt[Nmut_aa==1],ggplot2::aes(fitness,sigma)) + 
      ggplot2::geom_hex() + 
      ggplot2::scale_y_log10() + 
      ggplot2::coord_cartesian(ylim=c(0.01,1))
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_5_sigma_vs_fitness_singles.png"), d, width = 5, height = 5)
  }

  #### doubles
  #uncorrected fitness
  fitness_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]_uncorr"),colnames(doubles_dt))]
  sigma_rx <- sqrt(doubles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_uncorr"),colnames(doubles_dt))]^2 + 
                    matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
  # sigma_s2 <- random_effect_model(fitness_rx,sigma_rx)
  # doubles_dt[,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2 + sigma_s2),na.rm=T)/rowSums(1/(sigma_rx^2 + sigma_s2),na.rm=T)]
  doubles_dt[,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  # doubles_dt[,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2+sigma_s2),na.rm=T))]
  doubles_dt[,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  if(report){
    d <- ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_uncorr,sigma_uncorr)) + 
      ggplot2::geom_hex() + 
      ggplot2::scale_y_log10() + 
      ggplot2::coord_cartesian(ylim=c(0.01,1))
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_5_sigma_vs_fitness_doubles_uncorr.png"), d, width = 5, height = 5)
  }

  #conditioned fitness
  if(bayesian_double_fitness){
    fitness_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]_cond"),colnames(doubles_dt))]
    sigma_rx <- sqrt(doubles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_cond"),colnames(doubles_dt))]^2 + 
                      matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
    # sigma_s2 <- random_effect_model(fitness_rx,sigma_rx)
    # doubles_dt[,fitness_cond := rowSums(fitness_rx/(sigma_rx^2 + sigma_s2),na.rm=T)/rowSums(1/(sigma_rx^2 + sigma_s2),na.rm=T)]
    doubles_dt[,fitness_cond := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
    # doubles_dt[,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2+sigma_s2),na.rm=T))]
    doubles_dt[,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
    if(report){
      d <- ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_cond,sigma_cond)) + 
        ggplot2::geom_hex() + 
        ggplot2::scale_y_log10() + 
        ggplot2::coord_cartesian(ylim=c(0.01,1))
      ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_5_sigma_vs_fitness_doubles_cond.png"), d, width = 5, height = 5)
    }
  }

  #Plot to compare double fitness estimates
  if(report){
    p1<-ggplot2::ggplot(doubles_dt,ggplot2::aes(mean_count,fitness_uncorr)) + 
      ggplot2::geom_hex() + 
      ggplot2::scale_x_log10() +
      ggplot2::scale_fill_continuous(trans="log10")
    p3<-ggplot2::ggplot(doubles_dt[between(bin_count,2,8)],ggplot2::aes(fitness_uncorr,..scaled..,color=factor(bin_count))) +
      ggplot2::geom_density(adjust=1)
    p5<-ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_uncorr,sigma_uncorr)) + 
      ggplot2::geom_hex() + 
      ggplot2::scale_y_log10() +
      ggplot2::coord_cartesian(ylim = c(0.05,1.5))
    if(bayesian_double_fitness){
      p2<-ggplot2::ggplot(doubles_dt,ggplot2::aes(mean_count,fitness_cond)) + 
        ggplot2::geom_hex()+ 
        ggplot2::scale_x_log10() +
        ggplot2::scale_fill_continuous(trans="log10")
      p4<-ggplot2::ggplot(doubles_dt[between(bin_count,2,8)],ggplot2::aes(fitness_cond,..scaled..,color=factor(bin_count))) +
        ggplot2::geom_density(adjust=1)
      p6<-ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_cond,sigma_cond)) + 
        ggplot2::geom_hex()+ 
        ggplot2::scale_y_log10() +
        ggplot2::coord_cartesian(ylim = c(0.05,1.5))
    }
    ggplot2::theme_set(ggplot2::theme_minimal())
    #Plot
    if(bayesian_double_fitness){
      d <- cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5,p6),nrow=3)
      rm(p1,p2,p3,p4,p5,p6)
      ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_5_doubles_fitness_estimates.png"), d, width = 10, height = 10)
    }else{
      d <- cowplot::plot_grid(plotlist = list(p1,p3,p5),nrow=3)
      rm(p1,p3,p5)
      ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_5_doubles_fitness_estimates.png"), d, width = 5, height = 10)
    }
  }

  #Plot fitness values against each other
  if(report & bayesian_double_fitness){
    set.seed(1)
    d <- GGally::ggpairs(doubles_dt[sample(.N,1000),.(fitness_uncorr,fitness_cond)])
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_5_doubles_fitness_estimates_scattermatrix.png"), d, width = 10, height = 10)
  }

  #Plot sigma values against each other
  if(report & bayesian_double_fitness){
    d <- GGally::ggpairs(doubles_dt[,.(sigma_uncorr,sigma_cond)])
    ggplot2::ggsave(file.path(report_outpath, "dimsum_stage_fitness_report_5_doubles_sigma_estimates_scattermatrix.png"), d, width = 10, height = 10)
  }

  message("Done")

  ### Output replicate data files
  ###########################

  message("Saving fitness estimates...")

  # define which variants have enough reads
  wt_dt[,is.reads0 := TRUE]
  singles_dt[mean_count >= min_mean_input_read_count,is.reads0 := TRUE]
  doubles_dt[mean_count >= min_mean_input_read_count,is.reads0 := TRUE]

  #Reformat columns
  if(dimsum_meta[["sequenceType"]]=="coding"){
    #Rename columns
    names(wt_dt)[names(wt_dt)=="merge_seq"] <- "nt_seq"
    #Remove unnecessary columns
    unnecessary_cols <- c(
      "var_fitness",
      "avg_sigma",
      "isNA",
      "avg_sigma_fit",
      "counts_for_bins",
      "bin_count")
    input_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% unnecessary_cols]]
    wt_dt <- wt_dt[,.SD,,.SDcols = names(wt_dt)[!names(wt_dt) %in% unnecessary_cols]]
    singles_dt <- singles_dt[,.SD,,.SDcols = names(singles_dt)[!names(singles_dt) %in% unnecessary_cols]]
    doubles_dt <- doubles_dt[,.SD,,.SDcols = names(doubles_dt)[!names(doubles_dt) %in% unnecessary_cols]]
  }else{
    #Rename columns
    names(input_dt)[names(input_dt)=="merge_seq"] <- "nt_seq"
    names(wt_dt)[names(wt_dt)=="merge_seq"] <- "nt_seq"
    #Remove unnecessary columns
    unnecessary_cols <- c(
      "aa_seq",
      "Nmut_aa",
      "Nmut_codons",
      "STOP",
      "var_fitness",
      "avg_sigma",
      "isNA",
      "avg_sigma_fit",
      "merge_seq",
      "counts_for_bins",
      "bin_count")
    input_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% unnecessary_cols]]
    wt_dt <- wt_dt[,.SD,,.SDcols = names(wt_dt)[!names(wt_dt) %in% unnecessary_cols]]
    singles_dt <- singles_dt[,.SD,,.SDcols = names(singles_dt)[!names(singles_dt) %in% unnecessary_cols]]
    doubles_dt <- doubles_dt[,.SD,,.SDcols = names(doubles_dt)[!names(doubles_dt) %in% unnecessary_cols]]
    #Rename WT_AA columns to WT_nt
    names(input_dt)[grep("^WT_AA", names(input_dt))] <- gsub("WT_AA", "WT_nt", names(input_dt)[grep("^WT_AA", names(input_dt))])
    names(wt_dt)[grep("^WT_AA", names(wt_dt))] <- gsub("WT_AA", "WT_nt", names(wt_dt)[grep("^WT_AA", names(wt_dt))])
    names(singles_dt)[grep("^WT_AA", names(singles_dt))] <- gsub("WT_AA", "WT_nt", names(singles_dt)[grep("^WT_AA", names(singles_dt))])
    names(doubles_dt)[grep("^WT_AA", names(doubles_dt))] <- gsub("WT_AA", "WT_nt", names(doubles_dt)[grep("^WT_AA", names(doubles_dt))])
  }

  #Rename objects
  all_variants <- input_dt
  wildtype <- wt_dt
  doubles <- doubles_dt
  if(dimsum_meta[["sequenceType"]]=="coding"){
    silent <- singles_dt[Nmut_aa==0]
    singles <- singles_dt[Nmut_aa==1]
    save(all_variants, wildtype, silent, singles, doubles, file = file.path(fitness_outpath, paste0(dimsum_meta[["projectName"]], '_fitness_replicates.RData')))
  }else{
    singles <- singles_dt
    save(all_variants, wildtype, singles, doubles, file = file.path(fitness_outpath, paste0(dimsum_meta[["projectName"]], '_fitness_replicates.RData')))
  }

  ### Output plain text files
  ###########################

  ##### finalize data.tables
  if(dimsum_meta[["sequenceType"]]=="coding"){
    silent <- singles_dt[Nmut_aa==0,.(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,fitness,sigma)]
    singles <- singles_dt[Nmut_aa==1,.(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,fitness,sigma)]
  }else{
    singles <- singles_dt[,.(Pos,WT_nt,Mut,Nmut_nt,mean_count,is.reads0,fitness,sigma)]
  }

  #for doubles #add single mutant fitness/sigma values to double mutant table
  doubles_dt[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
  doubles_dt[,sigma1 := singles[Pos == Pos1 & Mut == Mut1,sigma],.(Pos1,Mut1)]
  doubles_dt[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
  doubles_dt[,sigma2 := singles[Pos == Pos2 & Mut == Mut2,sigma],.(Pos2,Mut2)]
  if(dimsum_meta[["sequenceType"]]=="coding"){
    retained_cols <- c("Pos1","Pos2","WT_AA1","WT_AA2","Mut1","Mut2","Nmut_nt","Nmut_aa","Nmut_codons","STOP","mean_count","is.reads0",
      "fitness1","sigma1","fitness2","sigma2",
      "fitness_uncorr","sigma_uncorr",
      "fitness_cond","sigma_cond")
    doubles <- doubles_dt[,.SD,,.SDcols = names(doubles_dt)[names(doubles_dt) %in% retained_cols]]

    #Exclude variants with STOP codons from downstream fitness analyses
    wildtype[,is.fitness := TRUE]
    silent[,is.fitness := TRUE]
    singles[,is.fitness := !STOP]
    doubles[,is.fitness := !STOP]

    #write data to files
    write.table(x = wildtype, file = file.path(fitness_outpath, "fitness_wildtype.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = silent, file = file.path(fitness_outpath, "fitness_silent.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = singles, file = file.path(fitness_outpath, "fitness_singles.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = doubles, file = file.path(fitness_outpath, "fitness_doubles.txt"),
                quote = F,row.names = F, col.names = T)
  }else{
    retained_cols <- c("Pos1","Pos2","WT_nt1","WT_nt2","Mut1","Mut2","Nmut_nt","mean_count","is.reads0",
      "fitness1","sigma1","fitness2","sigma2",
      "fitness_uncorr","sigma_uncorr",
      "fitness_cond","sigma_cond")
    doubles <- doubles_dt[,.SD,,.SDcols = names(doubles_dt)[names(doubles_dt) %in% retained_cols]]

    #Exclude variants with STOP codons from downstream fitness analyses
    wildtype[,is.fitness := TRUE]
    singles[,is.fitness := TRUE]
    doubles[,is.fitness := TRUE]

    #write data to files
    write.table(x = wildtype, file = file.path(fitness_outpath, "fitness_wildtype.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = singles, file = file.path(fitness_outpath, "fitness_singles.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = doubles, file = file.path(fitness_outpath, "fitness_doubles.txt"),
                quote = F,row.names = F, col.names = T)
  }

  message("Done")

}
