
#dimsum_stage_merge_report
#
# Generate final summary plots for all samples.
#
# dimsum_meta: an experiment metadata object (required)
# report_outpath: final report output path (required)
#
# Returns: an updated experiment metadata object
#
dimsum_stage_merge_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))
  #Final statistics
  merge_df <- dimsum_meta[['exp_design']]
  merge_df$pairname <- sapply(strsplit(merge_df$aligned_pair, '.split'), '[', 1)
  merge_df_collapse <- ddply(merge_df, "pairname", summarise, 
    total_read_pairs = sum(total_read_pairs), 
    usearch_merged = sum(usearch_merged), 
    cutadapt_pairs_too_short = sum(cutadapt_pairs_too_short),
    filter_incorrect_length = sum(filter_incorrect_length),
    too_many_aa_mutations = sum(too_many_aa_mutations))
  merge_df_collapse$usearch_not_merged <- merge_df_collapse$total_read_pairs-merge_df_collapse$cutadapt_pairs_too_short-merge_df_collapse$usearch_merged
  merge_df_collapse$retained <- merge_df_collapse$total_read_pairs-merge_df_collapse$cutadapt_pairs_too_short-merge_df_collapse$usearch_not_merged-merge_df_collapse$filter_incorrect_length-merge_df_collapse$too_many_aa_mutations
  merge_df_collapse_perc <- merge_df_collapse
  merge_df_collapse_perc[,4:8] <- merge_df_collapse_perc[,4:8]/merge_df_collapse_perc$total_read_pairs*100
  merge_df_collapse_perc <- merge_df_collapse_perc[,c(1,4:8)]
  merge_df_collapse <- merge_df_collapse[,c(1,4:8)]
  #Plot 1: final results for all samples (total reads)
  plot_df <- melt(merge_df_collapse, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Total read pairs (variants)", title = paste0("Read pairs to filtered variant statistics"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_variantcounts.png')), d, width=12, height=8)
  #Plot 2: final results for all samples (percentages)
  plot_df <- melt(merge_df_collapse_perc, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Percentage of initial read pairs (variants)", title = paste0("Read pairs to filtered variant statistics (percentage)"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_variantpercentages.png')), d, width=12, height=8)
  #Plot 3: AA mutation counts
  aa_mut_df <- data.frame(
    'AA_mut_0'=sapply(dimsum_meta[['aa_mutation_counts']], '[', '0'),
    'AA_mut_1'=sapply(dimsum_meta[['aa_mutation_counts']], '[', '1'),
    'AA_mut_2'=sapply(dimsum_meta[['aa_mutation_counts']], '[', '2'),
    'AA_mut_sum'=sapply(dimsum_meta[['aa_mutation_counts']], sum))
  aa_mut_df[is.na(aa_mut_df)] <- 0
  aa_mut_df$pairname <- sapply(strsplit(merge_df$aligned_pair, '.split'), '[', 1)
  aa_mut_df_collapse <- ddply(aa_mut_df, "pairname", summarise, 
    AA_mut_0 = sum(AA_mut_0), 
    AA_mut_1 = sum(AA_mut_1), 
    AA_mut_2 = sum(AA_mut_2),
    AA_mut_sum = sum(AA_mut_sum))
  aa_mut_df_collapse$AA_mut_3plus = aa_mut_df_collapse$AA_mut_sum-aa_mut_df_collapse$AA_mut_0-aa_mut_df_collapse$AA_mut_1-aa_mut_df_collapse$AA_mut_2
  aa_mut_df_collapse_perc = aa_mut_df_collapse
  aa_mut_df_collapse_perc[,c(2:5, 6)] <- aa_mut_df_collapse_perc[,c(2:5, 6)]/aa_mut_df_collapse_perc$AA_mut_sum*100
  aa_mut_df_collapse_perc <- aa_mut_df_collapse_perc[,c(1:4, 6)]
  aa_mut_df_collapse <- aa_mut_df_collapse[,c(1:4, 6)]
  #Plot
  plot_df <- melt(aa_mut_df_collapse, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Total variants", title = paste0("Variant amino acid mutation statistics"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_aamutationcounts.png')), d, width=12, height=8)
  #Plot 4: AA mutation percentages
  plot_df <- melt(aa_mut_df_collapse_perc, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Percentage of variants", title = paste0("Variant amino acid mutation statistics (percentage)"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_aamutationpercentages.png')), d, width=12, height=8)
  #Plot 5: Nucleotide mutation counts
  nuc_mut_df <- data.frame(
    'nuc_mut_0'=sapply(dimsum_meta[['nuc_mutation_counts']], '[', '0'),
    'nuc_mut_1'=sapply(dimsum_meta[['nuc_mutation_counts']], '[', '1'),
    'nuc_mut_2'=sapply(dimsum_meta[['nuc_mutation_counts']], '[', '2'),
    'nuc_mut_sum'=sapply(dimsum_meta[['nuc_mutation_counts']], sum))
  nuc_mut_df[is.na(nuc_mut_df)] <- 0
  nuc_mut_df$pairname <- sapply(strsplit(merge_df$aligned_pair, '.split'), '[', 1)
  nuc_mut_df_collapse <- ddply(nuc_mut_df, "pairname", summarise, 
    nuc_mut_0 = sum(nuc_mut_0), 
    nuc_mut_1 = sum(nuc_mut_1), 
    nuc_mut_2 = sum(nuc_mut_2),
    nuc_mut_sum = sum(nuc_mut_sum))
  nuc_mut_df_collapse$nuc_mut_3plus = nuc_mut_df_collapse$nuc_mut_sum-nuc_mut_df_collapse$nuc_mut_0-nuc_mut_df_collapse$nuc_mut_1-nuc_mut_df_collapse$nuc_mut_2
  nuc_mut_df_collapse_perc = nuc_mut_df_collapse
  nuc_mut_df_collapse_perc[,c(2:5, 6)] <- nuc_mut_df_collapse_perc[,c(2:5, 6)]/nuc_mut_df_collapse_perc$nuc_mut_sum*100
  nuc_mut_df_collapse_perc <- nuc_mut_df_collapse_perc[,c(1:4, 6)]
  nuc_mut_df_collapse <- nuc_mut_df_collapse[,c(1:4, 6)]
  #Plot
  plot_df <- melt(nuc_mut_df_collapse, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Total variants", title = paste0("Variant nucleotide mutation statistics"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_nucmutationcounts.png')), d, width=12, height=8)
  #Plot 6: Nucleotide mutation percentages
  plot_df <- melt(nuc_mut_df_collapse_perc, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Percentage of variants", title = paste0("Variant nucleotide mutation statistics (percentage)"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_nucmutationpercentages.png')), d, width=12, height=8)
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  return(dimsum_meta_new)
}

