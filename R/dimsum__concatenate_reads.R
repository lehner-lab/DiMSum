
#' dimsum__concatenate_reads
#'
#' Concatentate reads (without reverse complementing second read in pair).
#'
#' @param input_FASTQ1 Path to first read FASTQ file (required)
#' @param input_FASTQ2 Path to second read FASTQ file (required)
#' @param output_FASTQ Path to output FASTQ file (required)
#' @param output_REPORT Path to report file (required)
#' @param min_qual Minimum observed base quality to retain read pair (required)
#' @param max_ee Maximum number of expected errors to retain read pair (required)
#' @param min_len Discard pair if either read is shorter than this (required)
#'
#' @return Nothing
#' @export
dimsum__concatenate_reads <- function(
  input_FASTQ1,
  input_FASTQ2,
  output_FASTQ,
  output_REPORT,
  min_qual,
  max_ee,
  min_len
  ){
  #Alignment statistics
  a_stats <- list()
  a_stats[['Pairs']] <- 0
  a_stats[['Merged']] <- 0
  a_stats[['Alignments_zero_diffs']] <- 0
  a_stats[['Too_many_diffs']] <- 0
  a_stats[['Exp.errs._too_high']] <- 0
  a_stats[['Min_Q_too_low']] <- 0
  a_stats[['Fwd_too_short']] <- 0
  a_stats[['Rev_too_short']] <- 0
  a_stats[['merged_lengths']] <- list()
  #Input FASTQ files
  fq1 <- ShortRead::readFastq(input_FASTQ1)
  fq2 <- ShortRead::readFastq(input_FASTQ2)
  #Update statistics
  a_stats[['Pairs']] <- length(fq1)
  #Read lengths
  fq1_lengths <- IRanges::width(ShortRead::sread(fq1))
  fq2_lengths <- IRanges::width(ShortRead::sread(fq2))
  #Update statistics
  a_stats[['Fwd_too_short']] <- sum(fq1_lengths<min_len)
  a_stats[['Rev_too_short']] <- sum(fq2_lengths<min_len)
  #Subset to sequences at least 64 bp long
  fq1 <- fq1[fq1_lengths>=min_len & fq2_lengths>=min_len]
  fq2 <- fq2[fq1_lengths>=min_len & fq2_lengths>=min_len]
  #Update statistics
  a_stats[['Alignments_zero_diffs']] <- sum(fq1_lengths>=min_len & fq2_lengths>=min_len)
  #Read quality matrices
  qmat1 <- as(Biostrings::quality(fq1), "matrix")
  qmat2 <- as(Biostrings::quality(fq2), "matrix")
  #Number of bases with qualities less than minimum specified?
  non_merge_num_bases_too_low_qual = apply(qmat1<min_qual, 1, sum, na.rm = T) + apply(qmat2<min_qual, 1, sum, na.rm = T)
  #Update statistics
  a_stats[['Min_Q_too_low']] <- sum(non_merge_num_bases_too_low_qual!=0)
  #Subset to sequences with all base qualities not less than specified
  fq1 <- fq1[non_merge_num_bases_too_low_qual==0]
  fq2 <- fq2[non_merge_num_bases_too_low_qual==0]
  #Read error probability matrices
  emat1 <- 10^(qmat1[non_merge_num_bases_too_low_qual==0,]/(-10))
  emat2 <- 10^(qmat2[non_merge_num_bases_too_low_qual==0,]/(-10))
  #Expected number of read errors
  exp_num_read_errors <- apply(emat1, 1, sum, na.rm = T) + apply(emat2, 1, sum, na.rm = T)
  #Update statistics
  a_stats[['Exp.errs._too_high']] <- sum(exp_num_read_errors>max_ee)
  #Subset to sequences with less than specified expected number of read errors
  fq1 <- fq1[exp_num_read_errors<=max_ee]
  fq2 <- fq2[exp_num_read_errors<=max_ee]
  #Concatenate sequence
  fqc <- ShortRead::ShortReadQ(
    sread = Biostrings::DNAStringSet(paste0(as.character(ShortRead::sread(fq1)), as.character(ShortRead::sread(fq2)))), 
    quality = Biostrings::BStringSet(paste0(as.character(Biostrings::quality(Biostrings::quality(fq1))), as.character(Biostrings::quality(Biostrings::quality(fq2))))), 
    id = ShortRead::id(fq1))
  #Write to file
  ShortRead::writeFastq(fqc, file = output_FASTQ, compress = F)
  #Update statistics
  a_stats[['Merged']] = length(fqc)
  if(a_stats[['Merged']]!=0){
    a_stats[['Merged_length_min']] = min(IRanges::width(ShortRead::sread(fqc)))
    a_stats[['Merged_length_low']] = as.numeric(quantile(IRanges::width(ShortRead::sread(fqc)), 0.25))
    a_stats[['Merged_length_median']] = as.numeric(median(IRanges::width(ShortRead::sread(fqc))))
    a_stats[['Merged_length_high']] = as.numeric(quantile(IRanges::width(ShortRead::sread(fqc)), 0.75))
    a_stats[['Merged_length_max']] = max(IRanges::width(ShortRead::sread(fqc)))
  }else{
    a_stats[['Merged_length_min']] = NA
    a_stats[['Merged_length_low']] = NA
    a_stats[['Merged_length_median']] = NA
    a_stats[['Merged_length_high']] = NA
    a_stats[['Merged_length_max']] = NA
  }
  #Report
  report_list <- list()
  report_list <- append(report_list, paste0(
    '\nMerge\n\tFwd ', 
    input_FASTQ1, 
    '\n\tRev ', 
    input_FASTQ2, 
    '\n\tKeep read labels\n\t', 
    as.character(a_stats['Merged']), 
    ' / ', 
    as.character(a_stats['Pairs']), 
    ' pairs merged\n\n'))
  report_list <- append(report_list, 'Merged length distribution:\n')
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_min']], '  Min\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_low']], '  Low quartile\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_median']], '  Median\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_high']], '  High quartile\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_max']], '  Max\n\nTotals:\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Pairs']], '  Pairs ()\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged']], '  Merged (, )\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Alignments_zero_diffs']], '  Alignments with zero diffs ()\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Too_many_diffs']], '  Too many diffs (> 1) ()\n'))
  report_list <- append(report_list, paste0('\t ', 0, '  Fwd tails Q <= 2 trimmed ()\n'))
  report_list <- append(report_list, paste0('\t ', 0, '  Rev tails Q <= 2 trimmed ()\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Fwd_too_short']], '  Fwd too short (< ', min_len, ') after tail trimming ()\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Rev_too_short']], '  Rev too short (< ', min_len, ') after tail trimming ()\n'))
  report_list <- append(report_list, paste0('\t ', 0, '  No alignment found ()\n'))
  report_list <- append(report_list, paste0('\t ', 0, '  Alignment too short ( ) ()\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Exp.errs._too_high']], '  Exp.errs. too high (max=', max_ee, ') ()\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Min_Q_too_low']], '  Min Q too low (<', min_qual, ') ()\n'))
  write(paste0(unlist(report_list), collapse = ""), file = output_REPORT, sep = "")
}




