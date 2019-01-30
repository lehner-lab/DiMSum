
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
  #Final statistics
  merge_df <- dimsum_meta[['exp_design']]
  merge_df[,'pairname'] <- sapply(strsplit(merge_df[,'aligned_pair'], '.split'), '[', 1)
  #Plot 1: AA mutation counts
  aa_subst_df <- data.frame(
    'AA_subst_0'=sapply(dimsum_meta[['aa_subst_counts']], '[', '0'),
    'AA_subst_1'=sapply(dimsum_meta[['aa_subst_counts']], '[', '1'),
    'AA_subst_2'=sapply(dimsum_meta[['aa_subst_counts']], '[', '2'),
    'AA_subst_sum'=sapply(dimsum_meta[['aa_subst_counts']], sum),
    'AA_indel_sum'=sapply(dimsum_meta[['aa_indel_counts']], sum))
  aa_subst_df[is.na(aa_subst_df)] <- 0
  aa_subst_df[,'pairname'] <- sapply(strsplit(merge_df[,'aligned_pair'], '.split'), '[', 1)
  aa_subst_df_collapse <- plyr::ddply(aa_subst_df, "pairname", plyr::summarise, 
    AA_subst_0 = sum(AA_subst_0), 
    AA_subst_1 = sum(AA_subst_1), 
    AA_subst_2 = sum(AA_subst_2),
    AA_subst_sum = sum(AA_subst_sum),
    AA_indel_sum = sum(AA_indel_sum))
  aa_subst_df_collapse[,'AA_subst_3plus'] <- aa_subst_df_collapse[,'AA_subst_sum']-aa_subst_df_collapse[,'AA_subst_0']-aa_subst_df_collapse[,'AA_subst_1']-aa_subst_df_collapse[,'AA_subst_2']
  aa_subst_df_collapse_perc <- aa_subst_df_collapse
  aa_subst_df_collapse_perc[,colnames(aa_subst_df_collapse_perc)[-1]] <- aa_subst_df_collapse_perc[,colnames(aa_subst_df_collapse_perc)[-1]]/(aa_subst_df_collapse_perc[,'AA_subst_sum']+aa_subst_df_collapse_perc[,'AA_indel_sum'])*100
  aa_subst_df_collapse_perc <- aa_subst_df_collapse_perc[,-which(colnames(aa_subst_df_collapse_perc) == "AA_subst_sum")]
  temp_colnames <- c(
    "0 AA substitutions",
    "1 AA substitutions",
    "2 AA substitutions",
    "AA indels",
    "3+ AA substitutions")
  colnames(aa_subst_df_collapse_perc)[-1] <- temp_colnames
  aa_subst_df_collapse <- aa_subst_df_collapse[,-which(colnames(aa_subst_df_collapse) == "AA_subst_sum")]
  colnames(aa_subst_df_collapse)[-1] <- temp_colnames
  #Plot
  plot_df <- reshape2::melt(aa_subst_df_collapse, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=temp_colnames[order(temp_colnames)])
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Total reads with mutation", title = paste0("Read amino acid mutation statistics"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_aamutationcounts.png')), d, width=12, height=8)
  #Plot 2: AA mutation percentages
  plot_df <- reshape2::melt(aa_subst_df_collapse_perc, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=temp_colnames[order(temp_colnames)])
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Percentage of reads with mutation", title = paste0("Read amino acid mutation statistics (percentage)"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_aamutationpercentages.png')), d, width=12, height=8)
  #Plot 5: Nucleotide mutation counts
  nuc_subst_df <- data.frame(
    'nuc_subst_0'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '0'),
    'nuc_subst_1'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '1'),
    'nuc_subst_2'=sapply(dimsum_meta[['nuc_subst_counts']], '[', '2'),
    'nuc_subst_sum'=sapply(dimsum_meta[['nuc_subst_counts']], sum),
    'nuc_indel_sum'=sapply(dimsum_meta[['nuc_indel_counts']], sum))
  nuc_subst_df[is.na(nuc_subst_df)] <- 0
  nuc_subst_df[,'pairname'] <- sapply(strsplit(merge_df[,'aligned_pair'], '.split'), '[', 1)
  nuc_subst_df_collapse <- plyr::ddply(nuc_subst_df, "pairname", plyr::summarise, 
    nuc_subst_0 = sum(nuc_subst_0), 
    nuc_subst_1 = sum(nuc_subst_1), 
    nuc_subst_2 = sum(nuc_subst_2),
    nuc_subst_sum = sum(nuc_subst_sum),
    nuc_indel_sum = sum(nuc_indel_sum))
  nuc_subst_df_collapse[,'nuc_subst_3plus'] <- nuc_subst_df_collapse[,'nuc_subst_sum']-nuc_subst_df_collapse[,'nuc_subst_0']-nuc_subst_df_collapse[,'nuc_subst_1']-nuc_subst_df_collapse[,'nuc_subst_2']
  nuc_subst_df_collapse_perc <- nuc_subst_df_collapse
  nuc_subst_df_collapse_perc[,colnames(nuc_subst_df_collapse_perc)[-1]] <- nuc_subst_df_collapse_perc[,colnames(nuc_subst_df_collapse_perc)[-1]]/(nuc_subst_df_collapse_perc[,'nuc_subst_sum']+nuc_subst_df_collapse_perc[,'nuc_indel_sum'])*100
  nuc_subst_df_collapse_perc <- nuc_subst_df_collapse_perc[,-which(colnames(nuc_subst_df_collapse_perc) == "nuc_subst_sum")]
  temp_colnames <- c(
    "0 nuc. substitutions",
    "1 nuc. substitutions",
    "2 nuc. substitutions",
    "nuc. indels",
    "3+ nuc. substitutions")
  colnames(nuc_subst_df_collapse_perc)[-1] <- temp_colnames
  nuc_subst_df_collapse <- nuc_subst_df_collapse[,-which(colnames(nuc_subst_df_collapse) == "nuc_subst_sum")]
  colnames(nuc_subst_df_collapse)[-1] <- temp_colnames
  #Plot
  plot_df <- reshape2::melt(nuc_subst_df_collapse, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=temp_colnames[order(temp_colnames)])
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Total reads with mutation", title = paste0("Read nucleotide mutation statistics"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_nucmutationcounts.png')), d, width=12, height=8)
  #Plot 6: Nucleotide mutation percentages
  plot_df <- reshape2::melt(nuc_subst_df_collapse_perc, id="pairname")
  plot_df[,'Mutation_type'] <- factor(plot_df[,'variable'], levels=temp_colnames[order(temp_colnames)])
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pairname, value)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mutation_type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample names", y = "Percentage of reads with mutation", title = paste0("Read nucleotide mutation statistics (percentage)"))
  ggplot2::ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_nucmutationpercentages.png')), d, width=12, height=8)
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  return(dimsum_meta_new)
}

