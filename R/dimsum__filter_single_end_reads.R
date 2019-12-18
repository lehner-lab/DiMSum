
#' dimsum__filter_single_end_reads
#'
#' Concatentate reads (without reverse complementing second read in pair).
#'
#' @param input_FASTQ Path to first read FASTQ file (required)
#' @param output_FASTQ Path to output FASTQ file (required)
#' @param output_REPORT Path to report file (required)
#' @param min_qual Minimum observed base quality to retain read pair (required)
#' @param max_ee Maximum number of expected errors to retain read pair (required)
#' @param min_len Discard pair if either read is shorter than this (required)
#'
#' @return Nothing
#' @export
dimsum__filter_single_end_reads <- function(
  input_FASTQ,
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
  a_stats[['merged_lengths']] <- c()

  #Process FASTQ files
  initial_write <- TRUE #records written to output file already?
  yield_size <- 1e6
  f <- ShortRead::FastqStreamer(input_FASTQ, n=yield_size)
  #Read input FASTQ files in chunks
  while(length(fq <- ShortRead::yield(f))){
    #Update statistics
    a_stats[['Pairs']] <- a_stats[['Pairs']] + length(fq)
    #Read lengths
    fq_lengths <- IRanges::width(ShortRead::sread(fq))
    #Update statistics
    a_stats[['Fwd_too_short']] <- a_stats[['Fwd_too_short']] + sum(fq_lengths<min_len)
    #Subset to sequences at least 64 bp long
    fq <- fq[fq_lengths>=min_len]
    #Update statistics
    a_stats[['Alignments_zero_diffs']] <- a_stats[['Alignments_zero_diffs']] + sum(fq_lengths>=min_len)
    #Read quality matrices
    qmat <- as(Biostrings::quality(fq), "matrix")

    #Go to next iteration if empty matrix
    if(sum(fq_lengths>=min_len)==0){next}
    #Number of bases with qualities less than minimum specified?
    if(sum(fq_lengths>=min_len)==1){
      non_merge_num_bases_too_low_qual <- sum(qmat<min_qual, na.rm = T)
    }else{
      non_merge_num_bases_too_low_qual <- apply(qmat<min_qual, 1, sum, na.rm = T)
    }
    #Update statistics
    a_stats[['Min_Q_too_low']] <- a_stats[['Min_Q_too_low']] + sum(non_merge_num_bases_too_low_qual!=0)
    #Subset to sequences with all base qualities not less than specified
    fq <- fq[non_merge_num_bases_too_low_qual==0]
    #Read error probability matrices
    emat <- 10^(qmat[non_merge_num_bases_too_low_qual==0,]/(-10))

    #Go to next iteration if empty matrix
    if(sum(non_merge_num_bases_too_low_qual==0)==0){next}
    #Expected number of read errors
    if(sum(non_merge_num_bases_too_low_qual==0)==1){
      exp_num_read_errors <- sum(emat, na.rm = T)
    }else{
      exp_num_read_errors <- apply(emat, 1, sum, na.rm = T)
    }
    #Update statistics
    a_stats[['Exp.errs._too_high']] <- a_stats[['Exp.errs._too_high']] + sum(exp_num_read_errors>max_ee)
    #Subset to sequences with less than specified expected number of read errors
    fq <- fq[exp_num_read_errors<=max_ee]

    #Go to next iteration if no fastq records remain
    if(sum(exp_num_read_errors<=max_ee)==0){next}
    #Write to file
    dimsum__writeFastq(shortreads = fq, outputFile = output_FASTQ, initial_write = initial_write)
    initial_write <- FALSE
    #Update statistics
    a_stats[['Merged']] <- a_stats[['Merged']] + length(fq)
    a_stats[['merged_lengths']] <- c(a_stats[['merged_lengths']], IRanges::width(ShortRead::sread(fq)))
  }
  
  #Write empty FASTQ file if no records written thus far
  if(initial_write){
    dimsum__writeFastq(shortreads = ShortRead::ShortReadQ(), outputFile = output_FASTQ, initial_write = initial_write)
  }

  #Update length statistics
  if(a_stats[['Merged']]!=0){
    a_stats[['Merged_length_min']] <- min(a_stats[['merged_lengths']])
    a_stats[['Merged_length_low']] <- as.numeric(quantile(a_stats[['merged_lengths']], 0.25))
    a_stats[['Merged_length_median']] <- as.numeric(median(a_stats[['merged_lengths']]))
    a_stats[['Merged_length_high']] <- as.numeric(quantile(a_stats[['merged_lengths']], 0.75))
    a_stats[['Merged_length_max']] <- max(a_stats[['merged_lengths']])
  }else{
    a_stats[['Merged_length_min']] <- NA
    a_stats[['Merged_length_low']] <- NA
    a_stats[['Merged_length_median']] <- NA
    a_stats[['Merged_length_high']] <- NA
    a_stats[['Merged_length_max']] <- NA
  }

  #Report
  report_list <- list()
  report_list <- append(report_list, paste0(
    '\nMerge\n\tFwd ', 
    input_FASTQ, 
    '\n\tRev ', 
    input_FASTQ, 
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




