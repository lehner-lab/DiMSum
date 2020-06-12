
#' dimsum__vsearch_report
#'
#' Generate VSEARCH summary plots for all samples.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param report_outpath VSEARCH report output path (required)
#'
#' @return an updated experiment metadata object
#' @export
dimsum__vsearch_report <- function(
  dimsum_meta,
  report_outpath
  ){

  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))

  #Input files
  vsearch_files <- file.path(dimsum_meta[['exp_design']][,'aligned_pair_directory'], gsub('.vsearch$', '.report', dimsum_meta[['exp_design']][,'aligned_pair']))
  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = vsearch_files,
    execute = TRUE,
    exit = FALSE)

  #Get cutadapt results for all read pairs
  vsearch_list <- list()
  for(i in 1:length(vsearch_files)){
    temp_out <- readLines(vsearch_files[i])
    vsearch_list[[i]] <- list()
    vsearch_list[[i]][['vsearch_merge_length_low_quartile']] <- as.integer(unlist(strsplit(temp_out[grep('Low quartile', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_merge_length_median']] <- as.integer(unlist(strsplit(temp_out[grep('Median', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_merge_length_high_quartile']] <- as.integer(unlist(strsplit(temp_out[grep('High quartile', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_total_read_pairs']] <- as.integer(unlist(strsplit(temp_out[grep('Pairs', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_merged']] <- as.integer(unlist(strsplit(temp_out[grep('Merged$', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_too_short']] <- as.integer(unlist(strsplit(temp_out[grep('Too short', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_no_alignment_found']] <- as.integer(unlist(strsplit(temp_out[grep('No alignment found', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_too_many_diffs']] <- as.integer(unlist(strsplit(temp_out[grep('Too many diffs', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_overlap_too_short']] <- as.integer(unlist(strsplit(temp_out[grep('Overlap too short', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_exp_errs_too_high']] <- as.integer(unlist(strsplit(temp_out[grep('Exp.errs. too high', temp_out)], ' '))[2])
    vsearch_list[[i]][['vsearch_min_Q_too_low']] <- as.integer(unlist(strsplit(temp_out[grep('Min Q too low', temp_out)], ' '))[2])
    #Replace any missing values with zero
    for(j in 1:length(vsearch_list[[i]])){
      if(length(vsearch_list[[i]][[j]])==0){
        vsearch_list[[i]][[j]] <- 0
      }
    }
  }
  vsearch_df <- as.data.frame(apply(as.data.frame(do.call('rbind', vsearch_list)), 2, unlist))
  #Merge experimental design with VSEARCH report statistics
  vsearch_df <- cbind(dimsum_meta[['exp_design']][,c('aligned_pair', 'total_read_pairs')], vsearch_df)
  vsearch_df[,'cutadapt_not_written'] <- vsearch_df[,'total_read_pairs'] - vsearch_df[,'vsearch_total_read_pairs']
  #Plot 1: read pair count statistics
  vsearch_df[,'pairname'] <- dimsum__plot_samplename(sapply(strsplit(vsearch_df[,'aligned_pair'], '.split'), '[', 1))
  vsearch_df_collapse <- plyr::ddply(vsearch_df, "pairname", plyr::summarise, 
    total_read_pairs = sum(total_read_pairs), 
    vsearch_aligned = sum(vsearch_merged), 
    vsearch_too_short = sum(vsearch_too_short), 
    vsearch_no_alignment_found = sum(vsearch_no_alignment_found), 
    vsearch_too_many_diffs = sum(vsearch_too_many_diffs), 
    vsearch_overlap_too_short = sum(vsearch_overlap_too_short), 
    vsearch_exp_errs_too_high = sum(vsearch_exp_errs_too_high),
    vsearch_min_Q_too_low = sum(vsearch_min_Q_too_low),
    cutadapt_not_written = sum(cutadapt_not_written)
    )
  vsearch_df_collapse_perc <- vsearch_df_collapse
  vsearch_df_collapse_perc[,3:10] <- as.data.frame(t(scale(t(vsearch_df_collapse_perc[,3:10]), center = F, scale = vsearch_df_collapse_perc[,'total_read_pairs'])))*100
  vsearch_df_collapse_perc <- vsearch_df_collapse_perc[,c(1, 3:10)]
  #Plot
  plot_df <- reshape2::melt(vsearch_df_collapse_perc, id="pairname")
  plot_df[,'Alignment_status'] <- factor(plot_df[,'variable'])
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Alignment_status)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Read pairs (percentage)")#, title = paste0("Read pair alignment statistics"))
  dimsum__save_png(file.path(report_outpath, paste0('dimsum__vsearch_report_paircounts.png')), d, width=12, height=8)
  #Plot2: read pair merge length statistics
  plot_df <- reshape2::melt(vsearch_df[,grep('pairname|vsearch_merge_', colnames(vsearch_df))], id="pairname")
  plot_df[,'Length_quantile'] <- factor(plot_df[,'variable'])
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_boxplot(ggplot2::aes(color = Length_quantile)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    # ggplot2::coord_cartesian(ylim = c(0, max(plot_df[,'value'], na.rm = T))) +
    ggplot2::labs(x = "Sample names", y = "Aligned length (bp)")#, title = paste0("Aligned length distributions (all splits)"))
  dimsum__save_png(file.path(report_outpath, paste0('dimsum__vsearch_report_mergedlength.png')), d, width=12, height=6)

  #Render report
  dimsum__render_report(dimsum_meta = dimsum_meta)

  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']] <- cbind(dimsum_meta_new[['exp_design']], vsearch_df[,3:14])
  return(dimsum_meta_new)
}

