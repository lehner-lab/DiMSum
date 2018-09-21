
#dimsum_stage_usearch_report
#
# Generate USEARCH summary plots for all samples.
#
# dimsum_meta: an experiment metadata object (required)
# report_outpath: USEARCH report output path (required)
#
# Returns: an updated experiment metadata object
#
dimsum_stage_usearch_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))
  #Get cutadapt results for all read pairs
  usearch_files <- file.path(dimsum_meta[['exp_design']]$aligned_pair_directory, gsub('.usearch$', '.report', dimsum_meta[['exp_design']][,'aligned_pair']))
  usearch_list <- list()
  for(i in 1:length(usearch_files)){
    temp_out <- system(paste0("cat ", usearch_files[i]), intern=TRUE)
    usearch_list[[i]] <- list()
    usearch_list[[i]][['usearch_merge_length_low_quartile']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Low quartile', temp_out)], ' ')))[4])
    usearch_list[[i]][['usearch_merge_length_median']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Median', temp_out)], ' ')))[3])
    usearch_list[[i]][['usearch_merge_length_high_quartile']] <- as.integer(rev(unlist(strsplit(temp_out[grep('High quartile', temp_out)], ' ')))[4])
    usearch_list[[i]][['usearch_total_read_pairs']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Pairs', temp_out)], ' ')))[4])
    usearch_list[[i]][['usearch_merged']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Merged', temp_out)], ' ')))[5])
    usearch_list[[i]][['usearch_too_many_diffs']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Too many diffs', temp_out)], ' ')))[8])
    usearch_list[[i]][['usearch_fwd_too_short']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Fwd too short', temp_out)], ' ')))[11])
    usearch_list[[i]][['usearch_rev_too_short']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Rev too short', temp_out)], ' ')))[11])
    usearch_list[[i]][['usearch_no_alignment_found']] <- as.integer(rev(unlist(strsplit(temp_out[grep('No alignment found', temp_out)], ' ')))[6])
    usearch_list[[i]][['usearch_alignment_too_short']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Alignment too short', temp_out)], ' ')))[8])
    usearch_list[[i]][['usearch_exp_errs_too_high']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Exp.errs. too high', temp_out)], ' ')))[7])
    usearch_list[[i]][['usearch_min_Q_too_low']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Min Q too low', temp_out)], ' ')))[8])
    #Replace any missing values with zero
    for(j in 1:length(usearch_list[[i]])){
      if(length(usearch_list[[i]][[j]])==0){
        usearch_list[[i]][[j]] <- 0
      }
    }
  }
  usearch_df <- as.data.frame(apply(as.data.frame(do.call('rbind', usearch_list)), 2, unlist))
  #Merge experimental design with USEARCH report statistics
  usearch_df <- cbind(dimsum_meta[['exp_design']][,c('aligned_pair', 'total_read_pairs')], usearch_df)
  usearch_df$cutadapt_pairs_too_short <- usearch_df$total_read_pairs - usearch_df$usearch_total_read_pairs
  #Plot 1: read pair count statistics
  usearch_df$pairname <- sapply(strsplit(usearch_df$aligned_pair, '.split'), '[', 1)
  usearch_df_collapse <- ddply(usearch_df, "pairname", summarise, 
    total_read_pairs = sum(total_read_pairs), 
    usearch_merged = sum(usearch_merged), 
    usearch_too_many_diffs = sum(usearch_too_many_diffs), 
    usearch_fwd_too_short = sum(usearch_fwd_too_short), 
    usearch_rev_too_short = sum(usearch_rev_too_short), 
    usearch_no_alignment_found = sum(usearch_no_alignment_found), 
    usearch_alignment_too_short = sum(usearch_alignment_too_short), 
    usearch_exp_errs_too_high = sum(usearch_exp_errs_too_high),
    usearch_min_Q_too_low = sum(usearch_min_Q_too_low),
    cutadapt_pairs_too_short = sum(cutadapt_pairs_too_short)
    )
  usearch_df_collapse_perc <- usearch_df_collapse
  usearch_df_collapse_perc[,3:11] <- as.data.frame(t(scale(t(usearch_df_collapse_perc[,3:11]), center = F, scale = usearch_df_collapse_perc$total_read_pairs)))*100
  usearch_df_collapse_perc <- usearch_df_collapse_perc[,c(1, 3:11)]
  #Plot
  plot_df <- melt(usearch_df_collapse_perc, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Pairs (percentage)", title = paste0("Read pair alignment statistics"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_usearch_report_paircounts.png')), d, width=12, height=8)
  #Plot2: read pair merge length statistics
  plot_df <- melt(usearch_df[,grep('pairname|usearch_merge_', colnames(usearch_df))], id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_boxplot(aes(color = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_cartesian(ylim = c(0, max(plot_df$value))) +
    labs(x = "Sample names", y = "Merged length", title = paste0("Alignment length distributions (over splits)"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_usearch_report_mergedlength.png')), d, width=12, height=8)
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']] <- cbind(dimsum_meta_new[['exp_design']], usearch_df[,3:15])
  return(dimsum_meta_new)
}
