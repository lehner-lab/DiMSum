
#dimsum_stage_fastqc_report
#
# Generate FASTQC summary plots for all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# report_outpath: FASTQC report output path (required)
#
# Returns: an updated experiment metadata object
#
dimsum_stage_fastqc_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))
  #Get results for all fastq files
  for(col_name in c('pair1_fastqc', 'pair2_fastqc')){
    fastqc_files <- file.path(dimsum_meta[['exp_design']]$fastqc_directory, dimsum_meta[['exp_design']][,col_name])
    fastqc_list <- list()
    encoding <- ''
    for(f in fastqc_files){
      temp_out <- system(paste0("head -n ", 500, ' ', f), intern=TRUE)
      filename <- gsub('Filename\\t', '', temp_out[4])
      encoding <- gsub('Encoding\\t', '', temp_out[6])
      temp_nlines <- grep('>>END_MODULE', temp_out)[2]
      temp_out <- temp_out[c(13:(temp_nlines-1))]
      temp_out_data <- strsplit(temp_out[2:length(temp_out)], '\\t')
      temp_df <- as.data.frame(do.call('rbind', lapply(lapply(temp_out_data, '[', -1), as.numeric)))
      rownames(temp_df) <- sapply(temp_out_data, '[', 1)
      colnames(temp_df) <- unlist(strsplit(temp_out[1], '\\t'))[-1]
      fastqc_list[[filename]] <- temp_df
    }
    fastqc_df1 <- cbind.fill(lapply(fastqc_list, '[', 'Mean'))
    colnames(fastqc_df1) <- names(fastqc_list)
    fastqc_df1$base_position <- 1:length(rownames(fastqc_df1))
    fastqc_df2 <- cbind.fill(lapply(fastqc_list, '[', '10th Percentile'))
    colnames(fastqc_df2) <- names(fastqc_list)
    fastqc_df2$base_position <- 1:length(rownames(fastqc_df1))
    #Plot
    plot_df1 <- melt(fastqc_df1, id="base_position")
    plot_df1$statistic <- 'Mean'
    plot_df2 <- melt(fastqc_df2, id="base_position")
    plot_df2$statistic <- '10th Percentile'
    plot_df <- rbind(plot_df1, plot_df2)
    #Remove NAs
    plot_df <- plot_df[!is.na(plot_df$value),]
    d <- ggplot(plot_df, aes(base_position, value, color = variable)) +
      geom_line() +
      geom_hline(yintercept=c(20, 28), linetype = 2) +
      theme_bw() +
      coord_cartesian(ylim = c(0, max(plot_df$value))) +
      scale_x_continuous(
      breaks = (1:length(rownames(fastqc_df1)))[seq(1, length(rownames(fastqc_df1)), 5)],
      label = rownames(fastqc_df1)[seq(1, length(rownames(fastqc_df1)), 5)]) +
      labs(x = "Position in read (bp)", y = "Quality score", title = paste0("Quality scores across all bases (", encoding, ")"))
    d <- d + facet_wrap(~statistic, nrow=2, ncol=1)
    ggsave(file.path(report_outpath, paste0('dimsum_stage_fastqc_report_', col_name, '.png')), d, width=12, height=8)
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  return(dimsum_meta_new)
}

