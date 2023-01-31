
#' dimsum__filter_reads
#'
#' Concatenate reads (with or without reverse complementing second read in pair).
#'
#' @param input_FASTQ Path to input FASTQ file (required)
#' @param input_REPORT Path to input report file (required)
#' @param output_FASTQ Path to output FASTQ file (required)
#' @param output_REPORT Path to output report file (required)
#' @param min_qual Minimum observed base quality to retain read pair (required)
#'
#' @return Nothing
#' @export
dimsum__filter_reads <- function(
  input_FASTQ,
  input_REPORT,
  output_FASTQ,
  output_REPORT,
  min_qual
  ){
  #Alignment statistics
  a_stats <- list()
  a_stats[['Pairs']] <- 0 #from vsearch report
  a_stats[['Merged']] <- 0
  a_stats[['Too_short']] <- 0 #from vsearch report
  a_stats[['No_alignment_found']] <- 0 #from vsearch report
  a_stats[['Too_many_diffs']] <- 0 #from vsearch report
  a_stats[['Overlap_too_short']] <- 0 #from vsearch report
  a_stats[['Exp.errs._too_high']] <- 0 #from vsearch report
  a_stats[['Min_Q_too_low']] <- 0
  a_stats[['merged_lengths']] <- c() 

  #Get vsearch results
  temp_out <- readLines(input_REPORT)
  a_stats[['Pairs']] <- sum(as.integer(rev(unlist(strsplit(temp_out[grep('Pairs$', temp_out)], ' ')))[3]), na.rm = T)
  a_stats[['Too_short']] <- sum(as.integer(rev(unlist(strsplit(temp_out[grep('reads too short', temp_out)], ' ')))[7]), na.rm = T)
  too_few_kmers <- sum(as.integer(rev(unlist(strsplit(temp_out[grep('too few kmers found on same diagonal', temp_out)], ' ')))[9]), na.rm = T)
  multiple_potential <- sum(as.integer(rev(unlist(strsplit(temp_out[grep('multiple potential alignments', temp_out)], ' ')))[5]), na.rm = T)
  score_too_low <- sum(as.integer(rev(unlist(strsplit(temp_out[grep('alignment score too low', temp_out)], ' ')))[11]), na.rm = T)
  a_stats[['No_alignment_found']] <- sum(too_few_kmers, multiple_potential, score_too_low, na.rm = T)
  a_stats[['Too_many_diffs']] <- sum(as.integer(rev(unlist(strsplit(temp_out[grep('too many differences', temp_out)], ' ')))[5]), na.rm = T)
  a_stats[['Overlap_too_short']] <- sum(as.integer(rev(unlist(strsplit(temp_out[grep('overlap too short', temp_out)], ' ')))[5]), na.rm = T)
  a_stats[['Exp.errs._too_high']] <- sum(as.integer(rev(unlist(strsplit(temp_out[grep('expected error too high', temp_out)], ' ')))[6]), na.rm = T)

  #Process FASTQ files
  initial_write <- TRUE #records written to output file already?
  yield_size <- 1e6
  f1 <- ShortRead::FastqStreamer(input_FASTQ, n=yield_size)
  #Read input FASTQ files in chunks
  while(length(fq1 <- ShortRead::yield(f1))){
    #Read quality matrices
    qmat1 <- as(Biostrings::quality(fq1), "matrix")
    #Number of bases with qualities less than minimum specified?
    non_merge_num_bases_too_low_qual <- apply(qmat1<min_qual, 1, sum, na.rm = T)
    #Update statistics
    a_stats[['Min_Q_too_low']] <- a_stats[['Min_Q_too_low']] + sum(non_merge_num_bases_too_low_qual!=0)
    #Subset to sequences with all base qualities not less than specified
    fq1 <- fq1[non_merge_num_bases_too_low_qual==0]
    #Write to file
    dimsum__writeFastq(shortreads = fq1, outputFile = output_FASTQ, initial_write = initial_write)
    initial_write <- FALSE
    #Update statistics
    a_stats[['Merged']] <- a_stats[['Merged']] + length(fq1)
    a_stats[['merged_lengths']] <- c(a_stats[['merged_lengths']], IRanges::width(ShortRead::sread(fq1)))
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
  report_list <- append(report_list, 'Merged length distribution:\n')
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_min']], '  Min\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_low']], '  Low quartile\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_median']], '  Median\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_high']], '  High quartile\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged_length_max']], '  Max\n\nTotals:\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Pairs']], '  Pairs\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Merged']], '  Merged\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Too_short']], '  Too short\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['No_alignment_found']], '  No alignment found\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Too_many_diffs']], '  Too many diffs\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Overlap_too_short']], '  Overlap too short\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Exp.errs._too_high']], '  Exp.errs. too high\n'))
  report_list <- append(report_list, paste0('\t ', a_stats[['Min_Q_too_low']], '  Min Q too low\n'))
  write(paste0(unlist(report_list), collapse = ""), file = output_REPORT, sep = "")

  #Delete input FASTQ file
  suppressWarnings(temp_out <- file.remove(input_FASTQ))
}




