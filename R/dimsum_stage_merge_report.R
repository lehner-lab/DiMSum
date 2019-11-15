
#' dimsum_stage_merge_report
#'
#' Generate final summary plots for all samples.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param report_outpath final report output path (required)
#'
#' @return a single data.frame where empty rows are filled with NAs
#' @export
dimsum_stage_merge_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))

  #Variant processing statistics - nucleotide substitutions
  merge_df <- dimsum_meta[['exp_design']]
  merge_df[,'pairname'] <- sapply(strsplit(merge_df[,'aligned_pair'], '_t'), '[', 1)  
  nuc_subst_df <- data.frame(
    'nuc_subst_0'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '0'),
    'nuc_subst_1'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '1'),
    'nuc_subst_2'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '2'),
    'nuc_subst_sum'=sapply(dimsum_meta[['nuc_subst_counts']], sum),
    'nuc_tmsub_sum'=sapply(dimsum_meta[['nuc_tmsub_counts']], sum),
    'nuc_frbdn_sum'=sapply(dimsum_meta[['nuc_frbdn_counts']], sum),
    'nuc_const_sum'=sapply(dimsum_meta[['nuc_const_counts']], sum),    
    'nuc_indel_sum'=sapply(dimsum_meta[['nuc_indel_counts']], sum))
  nuc_subst_df[is.na(nuc_subst_df)] <- 0
  nuc_subst_df[,'pairname'] <- unique(sapply(strsplit(merge_df[,'aligned_pair'], '_t'), '[', 1))
  nuc_subst_df_collapse <- plyr::ddply(nuc_subst_df, "pairname", plyr::summarise, 
    nuc_subst_0 = sum(nuc_subst_0), 
    nuc_subst_1 = sum(nuc_subst_1), 
    nuc_subst_2 = sum(nuc_subst_2),
    nuc_subst_sum = sum(nuc_subst_sum),
    nuc_tmsub_sum = sum(nuc_tmsub_sum),
    nuc_frbdn_sum = sum(nuc_frbdn_sum),
    nuc_const_sum = sum(nuc_const_sum),
    nuc_indel_sum = sum(nuc_indel_sum))
  nuc_subst_df_collapse[,'nuc_subst_3plus'] <- nuc_subst_df_collapse[,'nuc_subst_sum']-nuc_subst_df_collapse[,'nuc_subst_0']-nuc_subst_df_collapse[,'nuc_subst_1']-nuc_subst_df_collapse[,'nuc_subst_2']
  nuc_subst_df_collapse_perc <- nuc_subst_df_collapse
  nuc_subst_df_collapse_perc[,colnames(nuc_subst_df_collapse_perc)[-1]] <- nuc_subst_df_collapse_perc[,colnames(nuc_subst_df_collapse_perc)[-1]]/(
    nuc_subst_df_collapse_perc[,'nuc_subst_sum'] + 
    nuc_subst_df_collapse_perc[,'nuc_tmsub_sum'] + 
    nuc_subst_df_collapse_perc[,'nuc_frbdn_sum'] + 
    nuc_subst_df_collapse_perc[,'nuc_const_sum'] + 
    nuc_subst_df_collapse_perc[,'nuc_indel_sum']) * 100
  nuc_subst_df_collapse_perc <- nuc_subst_df_collapse_perc[,-which(colnames(nuc_subst_df_collapse_perc) == "nuc_subst_sum")]
  temp_colnames <- c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "too many",
    "not permitted",
    "internal constant region",
    "indel",
    "3+ hamming dist.")
  colnames(nuc_subst_df_collapse_perc)[-1] <- temp_colnames
  nuc_subst_df_collapse <- nuc_subst_df_collapse[,-which(colnames(nuc_subst_df_collapse) == "nuc_subst_sum")]
  colnames(nuc_subst_df_collapse)[-1] <- temp_colnames
  
  #Plot 1: Nucleotide mutation counts
  plot_df <- reshape2::melt(nuc_subst_df_collapse, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "3+ hamming dist.",
    "too many",
    "not permitted",
    "internal constant region",
    "indel"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Total reads with mutation", title = paste0("Read nucleotide mutation statistics"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_nucmutationcounts.png')), d, width=12, height=8)
  
  #Plot 2: Nucleotide mutation percentages
  plot_df <- reshape2::melt(nuc_subst_df_collapse_perc, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "3+ hamming dist.",
    "too many",
    "not permitted",
    "internal constant region",
    "indel"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Percentage of reads with mutation", title = paste0("Read nucleotide mutation statistics (percentage)"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_nucmutationpercentages.png')), d, width=12, height=8)

  #Variant processing statistics - amino acid substitutions
  merge_df <- dimsum_meta[['exp_design']]
  merge_df[,'pairname'] <- sapply(strsplit(merge_df[,'aligned_pair'], '_t'), '[', 1)  
  aa_subst_df <- data.frame(
    'aa_subst_0'=sapply(dimsum_meta[['aa_subst_counts']], '[', '0'),
    'aa_subst_1'=sapply(dimsum_meta[['aa_subst_counts']], '[', '1'),
    'aa_subst_2'=sapply(dimsum_meta[['aa_subst_counts']], '[', '2'),
    'aa_subst_sum'=sapply(dimsum_meta[['aa_subst_counts']], sum),
    'nuc_tmsub_sum'=sapply(dimsum_meta[['nuc_tmsub_counts']], sum),
    'nuc_frbdn_sum'=sapply(dimsum_meta[['nuc_frbdn_counts']], sum),
    'nuc_const_sum'=sapply(dimsum_meta[['nuc_const_counts']], sum),    
    'nuc_indel_sum'=sapply(dimsum_meta[['nuc_indel_counts']], sum))
  aa_subst_df[is.na(aa_subst_df)] <- 0
  aa_subst_df[,'pairname'] <- unique(sapply(strsplit(merge_df[,'aligned_pair'], '_t'), '[', 1))
  aa_subst_df_collapse <- plyr::ddply(aa_subst_df, "pairname", plyr::summarise, 
    aa_subst_0 = sum(aa_subst_0), 
    aa_subst_1 = sum(aa_subst_1), 
    aa_subst_2 = sum(aa_subst_2),
    aa_subst_sum = sum(aa_subst_sum),
    nuc_tmsub_sum = sum(nuc_tmsub_sum),
    nuc_frbdn_sum = sum(nuc_frbdn_sum),
    nuc_const_sum = sum(nuc_const_sum),
    nuc_indel_sum = sum(nuc_indel_sum))
  aa_subst_df_collapse[,'aa_subst_3plus'] <- aa_subst_df_collapse[,'aa_subst_sum']-aa_subst_df_collapse[,'aa_subst_0']-aa_subst_df_collapse[,'aa_subst_1']-aa_subst_df_collapse[,'aa_subst_2']
  aa_subst_df_collapse_perc <- aa_subst_df_collapse
  aa_subst_df_collapse_perc[,colnames(aa_subst_df_collapse_perc)[-1]] <- aa_subst_df_collapse_perc[,colnames(aa_subst_df_collapse_perc)[-1]]/(
    aa_subst_df_collapse_perc[,'aa_subst_sum'] + 
    aa_subst_df_collapse_perc[,'nuc_tmsub_sum'] + 
    aa_subst_df_collapse_perc[,'nuc_frbdn_sum'] + 
    aa_subst_df_collapse_perc[,'nuc_const_sum'] + 
    aa_subst_df_collapse_perc[,'nuc_indel_sum']) * 100
  aa_subst_df_collapse_perc <- aa_subst_df_collapse_perc[,-which(colnames(aa_subst_df_collapse_perc) == "aa_subst_sum")]
  temp_colnames <- c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "too many",
    "not permitted",
    "internal constant region",
    "indel",
    "3+ hamming dist.")
  colnames(aa_subst_df_collapse_perc)[-1] <- temp_colnames
  aa_subst_df_collapse <- aa_subst_df_collapse[,-which(colnames(aa_subst_df_collapse) == "aa_subst_sum")]
  colnames(aa_subst_df_collapse)[-1] <- temp_colnames

  #Plot 3: Nucleotide mutation counts
  plot_df <- reshape2::melt(aa_subst_df_collapse, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "3+ hamming dist.",
    "too many",
    "not permitted",
    "internal constant region",
    "indel"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Total reads with mutation", title = paste0("Read amino acid mutation statistics"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_aamutationcounts.png')), d, width=12, height=8)
  
  #Plot 4: Nucleotide mutation percentages
  plot_df <- reshape2::melt(aa_subst_df_collapse_perc, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=c(
    "0 hamming dist.",
    "1 hamming dist.",
    "2 hamming dist.",
    "3+ hamming dist.",
    "too many",
    "not permitted",
    "internal constant region",
    "indel"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Percentage of reads with mutation", title = paste0("Read amino acid mutation statistics (percentage)"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_aamutationpercentages.png')), d, width=12, height=8)

  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  return(dimsum_meta_new)
}

