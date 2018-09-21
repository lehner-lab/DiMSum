
#dimsum_stage_cutadapt_report
#
# Generate cutadapt summary plots for all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# report_outpath: cutadapt report output path (required)
#
# Returns: an updated experiment metadata object
#
dimsum_stage_cutadapt_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))
  #Get cutadapt results for all read pairs
  cutadapt_files <- file.path(dimsum_meta[['exp_design']]$pair_directory, paste0(dimsum_meta[['exp_design']][,'pair1'], '.stdout'))
  cutadapt_read1_list <- list()
  cutadapt_read2_list <- list()
  total_reads_list <- list()
  for(i in 1:length(cutadapt_files)){
    temp_out <- system(paste0("cat ", cutadapt_files[i]), intern=TRUE)
    name_read1 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[2])
    name_read2 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
    temp_out <- temp_out[c(9:13, grep("Sequence: ", temp_out))]
    total_reads <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[1], ' ')))[1]))
    total_reads_list[[i]] <- total_reads
    total_read1_trimmed <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[2], ' ')))[2]))
    total_read2_trimmed <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[3], ' ')))[2]))
    total_read1_a3 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[6], ' ')))[2]))
    total_read1_a5 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[7], ' ')))[2]))
    total_read1_both <- total_read1_a3+total_read1_a5-total_read1_trimmed
    total_read2_a3 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[8], ' ')))[2]))
    total_read2_a5 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[9], ' ')))[2]))
    total_read2_both <- total_read2_a3+total_read2_a5-total_read2_trimmed
    #If linked adapter supplied
    if(grepl('Type: linked', temp_out[6]) & grepl('Type: linked', temp_out[7])){
      total_read1_a3 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[6], ' ')))[2]))
      total_read1_a5 <- total_read1_a3
      total_read1_both <- total_read1_a3
      total_read2_a3 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[7], ' ')))[2]))
      total_read2_a5 <- total_read2_a3
      total_read2_both <- total_read2_a3
    }
    cutadapt_read1_list[[name_read1]] <- c(
      total_read1_a5-total_read1_both, total_read1_a3-total_read1_both, total_read1_both, total_reads)
    cutadapt_read2_list[[name_read2]] <- c(
      total_read2_a5-total_read2_both, total_read2_a3-total_read2_both, total_read2_both, total_reads)
  }
  #First read
  cutadapt_read1_df <- as.data.frame(do.call('rbind', cutadapt_read1_list))
  colnames(cutadapt_read1_df) <- c('adapter5prime', 'adapter3prime', 'both', 'total_reads')
  cutadapt_read1_df$fastq <- sapply(strsplit(rownames(cutadapt_read1_df), '.split'), '[', 1)
  cutadapt_read1_df_collapse <- ddply(cutadapt_read1_df, "fastq", summarise, 
    adapter5prime = sum(adapter5prime), 
    adapter3prime = sum(adapter3prime), 
    both = sum(both), 
    total_reads = sum(total_reads))
  cutadapt_read1_df_collapse_perc <- cutadapt_read1_df_collapse
  cutadapt_read1_df_collapse_perc[,2:4] <- as.data.frame(t(scale(t(cutadapt_read1_df_collapse_perc[,2:4]), center = F, scale = cutadapt_read1_df_collapse_perc$total_reads)))*100
  cutadapt_read1_df_collapse_perc <- cutadapt_read1_df_collapse_perc[,1:4]
  #Plot
  plot_df <- melt(cutadapt_read1_df_collapse_perc, id="fastq")
  d <- ggplot(plot_df, aes(fastq, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "FASTQ files", y = "Reads with adapters (percentage)", title = paste0("Percentage first read trimmed with cutadapt"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_cutadapt_report_pair1.png')), d, width=12, height=8)
  #Second read
  cutadapt_read2_df <- as.data.frame(do.call('rbind', cutadapt_read2_list))
  colnames(cutadapt_read2_df) <- c('adapter5prime', 'adapter3prime', 'both', 'total_reads')
  cutadapt_read2_df$fastq <- sapply(strsplit(rownames(cutadapt_read2_df), '.split'), '[', 1)
  cutadapt_read2_df_collapse <- ddply(cutadapt_read2_df, "fastq", summarise, 
    adapter5prime = sum(adapter5prime), 
    adapter3prime = sum(adapter3prime), 
    both = sum(both), 
    total_reads = sum(total_reads))
  cutadapt_read2_df_collapse_perc <- cutadapt_read2_df_collapse
  cutadapt_read2_df_collapse_perc[,2:4] <- as.data.frame(t(scale(t(cutadapt_read2_df_collapse_perc[,2:4]), center = F, scale = cutadapt_read2_df_collapse_perc$total_reads)))*100
  cutadapt_read2_df_collapse_perc <- cutadapt_read2_df_collapse_perc[,1:4]
  #Plot
  plot_df <- melt(cutadapt_read2_df_collapse_perc, id="fastq")
  d <- ggplot(plot_df, aes(fastq, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "FASTQ files", y = "Reads with adapters (percentage)", title = paste0("Percentage second read trimmed with cutadapt"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_cutadapt_report_pair2.png')), d, width=12, height=8)
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']]$total_read_pairs <- as.numeric(unlist(total_reads_list))
  return(dimsum_meta_new)
}

