
#' dimsum__cutadapt_report
#'
#' Generate cutadapt summary plots for all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param report_outpath cutadapt report output path (required)
#'
#' @return an updated experiment metadata object
#' @export
dimsum__cutadapt_report <- function(
  dimsum_meta,
  report_outpath
  ){

  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))

  #Input files
  cutadapt_files <- file.path(dimsum_meta[['exp_design']][,'pair_directory'], paste0(dimsum_meta[['exp_design']][,'pair1'], '.stdout'))
  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = cutadapt_files,
    execute = TRUE,
    exit = FALSE)

  #Get cutadapt results for all read pairs
  cutadapt_read1_list <- list()
  cutadapt_read2_list <- list()
  total_reads_list <- list()
  for(i in 1:length(cutadapt_files)){
    if(dimsum_meta[["paired"]]){
      trim_list <- dimsum__parse_cutadapt_output(
        file_path = cutadapt_files[i], 
        ran_cutadapt = dimsum_meta[['exp_design']][i,'run_cutadapt'],
        ran_cutadapt_cutonly = dimsum_meta[['exp_design']][i,'run_cutadapt_cutonly'])
      total_reads_list[[i]] <- trim_list[['total_reads']]
      cutadapt_read1_list[[trim_list[['name_read1']]]] <- c(
        trim_list[['total_read1_a5']]-trim_list[['total_read1_both']], 
        trim_list[['total_read1_a3']]-trim_list[['total_read1_both']], 
        trim_list[['total_read1_both']], 
        trim_list[['total_reads']])
      cutadapt_read2_list[[trim_list[['name_read2']]]] <- c(
        trim_list[['total_read2_a5']]-trim_list[['total_read2_both']], 
        trim_list[['total_read2_a3']]-trim_list[['total_read2_both']], 
        trim_list[['total_read2_both']], 
        trim_list[['total_reads']])
    }else{
      trim_list <- dimsum__parse_cutadapt_output_single_end(
        file_path = cutadapt_files[i], 
        ran_cutadapt = dimsum_meta[['exp_design']][i,'run_cutadapt'],
        ran_cutadapt_cutonly = dimsum_meta[['exp_design']][i,'run_cutadapt_cutonly'])
      total_reads_list[[i]] <- trim_list[['total_reads']]
      cutadapt_read1_list[[trim_list[['name_read1']]]] <- c(
        trim_list[['total_read1_a5']]-trim_list[['total_read1_both']], 
        trim_list[['total_read1_a3']]-trim_list[['total_read1_both']], 
        trim_list[['total_read1_both']], 
        trim_list[['total_reads']])
    }
  }
  #First read
  cutadapt_read1_df <- as.data.frame(do.call('rbind', cutadapt_read1_list))
  colnames(cutadapt_read1_df) <- c('five_prime', 'three_prime', 'both', 'total_reads')
  cutadapt_read1_df[,'fastq'] <- sapply(strsplit(rownames(cutadapt_read1_df), '.split'), '[', 1)
  cutadapt_read1_df_collapse <- plyr::ddply(cutadapt_read1_df, "fastq", plyr::summarise, 
    five_prime = sum(five_prime), 
    three_prime = sum(three_prime), 
    both = sum(both), 
    total_reads = sum(total_reads))
  cutadapt_read1_df_collapse_perc <- cutadapt_read1_df_collapse
  cutadapt_read1_df_collapse_perc[,2:4] <- as.data.frame(t(scale(t(cutadapt_read1_df_collapse_perc[,2:4]), center = F, scale = cutadapt_read1_df_collapse_perc[,'total_reads'])))*100
  cutadapt_read1_df_collapse_perc <- cutadapt_read1_df_collapse_perc[,1:4]
  #Plot
  plot_df <- reshape2::melt(cutadapt_read1_df_collapse_perc, id="fastq")
  plot_df[,'Region_trimmed'] <- factor(plot_df[,'variable'])
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(fastq, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Region_trimmed)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "FASTQ files", y = "Reads trimmed (percentage)")#, title = paste0("Read 1 percentage constant region identified and trimmed"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum__cutadapt_report_pair1.png')), d, width=12, height=8)
  #Second read (if paired design)
  if(dimsum_meta[["paired"]]){
    cutadapt_read2_df <- as.data.frame(do.call('rbind', cutadapt_read2_list))
    colnames(cutadapt_read2_df) <- c('five_prime', 'three_prime', 'both', 'total_reads')
    cutadapt_read2_df[,'fastq'] <- sapply(strsplit(rownames(cutadapt_read2_df), '.split'), '[', 1)
    cutadapt_read2_df_collapse <- plyr::ddply(cutadapt_read2_df, "fastq", plyr::summarise, 
      five_prime = sum(five_prime), 
      three_prime = sum(three_prime), 
      both = sum(both), 
      total_reads = sum(total_reads))
    cutadapt_read2_df_collapse_perc <- cutadapt_read2_df_collapse
    cutadapt_read2_df_collapse_perc[,2:4] <- as.data.frame(t(scale(t(cutadapt_read2_df_collapse_perc[,2:4]), center = F, scale = cutadapt_read2_df_collapse_perc[,'total_reads'])))*100
    cutadapt_read2_df_collapse_perc <- cutadapt_read2_df_collapse_perc[,1:4]
    #Plot
    plot_df <- reshape2::melt(cutadapt_read2_df_collapse_perc, id="fastq")
    plot_df[,'Region_trimmed'] <- factor(plot_df[,'variable'])
    d <- ggplot2::ggplot(plot_df, ggplot2::aes(fastq, value)) +
      ggplot2::geom_col(ggplot2::aes(fill = Region_trimmed)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(x = "FASTQ files", y = "Reads trimmed (percentage)")#, title = paste0("Read 2 percentage constant region identified and trimmed"))
    ggplot2::ggsave(file.path(report_outpath, paste0('dimsum__cutadapt_report_pair2.png')), d, width=12, height=8)
  }

  #Render report
  dimsum__render_report(dimsum_meta = dimsum_meta)

  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']][,'total_read_pairs'] <- as.numeric(unlist(total_reads_list))
  return(dimsum_meta_new)
}

