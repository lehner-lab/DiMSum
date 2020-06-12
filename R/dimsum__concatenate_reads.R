
#' dimsum__concatenate_reads
#'
#' Concatenate reads (with or without reverse complementing second read in pair).
#'
#' @param input_FASTQ1 Path to first read FASTQ file (required)
#' @param input_FASTQ2 Path to second read FASTQ file (required)
#' @param output_FASTQ Path to output FASTQ file (required)
#' @param output_REPORT Path to report file (required)
#' @param min_qual Minimum observed base quality to retain read pair (required)
#' @param max_ee Maximum number of expected errors to retain read pair (required)
#' @param min_len Discard pair if either read is shorter than this (required)
#' @param reverse_complement_second_read Reverse complement second read before concatenation (default:FALSE)
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
  min_len,
  reverse_complement_second_read = FALSE
  ){
  #Alignment statistics
  a_stats <- list()
  a_stats[['Pairs']] <- 0
  a_stats[['Merged']] <- 0
  a_stats[['Too_short']] <- 0
  a_stats[['No_alignment_found']] <- 0 #NA
  a_stats[['Too_many_diffs']] <- 0 #NA
  a_stats[['Overlap_too_short']] <- 0 #NA
  a_stats[['Exp.errs._too_high']] <- 0
  a_stats[['Min_Q_too_low']] <- 0
  a_stats[['merged_lengths']] <- c()

  #Process FASTQ files
  initial_write <- TRUE #records written to output file already?
  yield_size <- 1e6
  f1 <- ShortRead::FastqStreamer(input_FASTQ1, n=yield_size)
  f2 <- ShortRead::FastqStreamer(input_FASTQ2, n=yield_size)
  #Read input FASTQ files in chunks
  while(length(fq1 <- ShortRead::yield(f1))){
    fq2 <- ShortRead::yield(f2)
    #Update statistics
    a_stats[['Pairs']] <- a_stats[['Pairs']] + length(fq1)
    #Read lengths
    fq1_lengths <- IRanges::width(ShortRead::sread(fq1))
    fq2_lengths <- IRanges::width(ShortRead::sread(fq2))
    #Update statistics
    a_stats[['Too_short']] <- a_stats[['Too_short']] + sum(fq1_lengths<min_len | fq2_lengths<min_len)
    #Subset to sequences at least 64 bp long
    fq1 <- fq1[fq1_lengths>=min_len & fq2_lengths>=min_len]
    fq2 <- fq2[fq1_lengths>=min_len & fq2_lengths>=min_len]

    #Read quality matrices
    qmat1 <- as(Biostrings::quality(fq1), "matrix")
    qmat2 <- as(Biostrings::quality(fq2), "matrix")
    #Number of bases with qualities less than minimum specified?
    non_merge_num_bases_too_low_qual <- apply(qmat1<min_qual, 1, sum, na.rm = T) + apply(qmat2<min_qual, 1, sum, na.rm = T)
    #Update statistics
    a_stats[['Min_Q_too_low']] <- a_stats[['Min_Q_too_low']] + sum(non_merge_num_bases_too_low_qual!=0)
    #Subset to sequences with all base qualities not less than specified
    fq1 <- fq1[non_merge_num_bases_too_low_qual==0]
    fq2 <- fq2[non_merge_num_bases_too_low_qual==0]
    #Read error probability matrices
    emat1 <- 10^(qmat1[non_merge_num_bases_too_low_qual==0,]/(-10))
    emat2 <- 10^(qmat2[non_merge_num_bases_too_low_qual==0,]/(-10))
    #Expected number of read errors
    exp_num_read_errors <- apply(emat1, 1, sum, na.rm = T) + apply(emat2, 1, sum, na.rm = T)
    #Update statistics
    a_stats[['Exp.errs._too_high']] <- a_stats[['Exp.errs._too_high']] + sum(exp_num_read_errors>max_ee)
    #Subset to sequences with less than specified expected number of read errors
    fq1 <- fq1[exp_num_read_errors<=max_ee]
    fq2 <- fq2[exp_num_read_errors<=max_ee]
    #Reverse complement second read if necessary
    if(reverse_complement_second_read){fq2 <- Biostrings::reverseComplement(fq2)}
    #Concatenate sequence
    fqc <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(paste0(as.character(ShortRead::sread(fq1)), as.character(ShortRead::sread(fq2)))), 
      quality = Biostrings::BStringSet(paste0(as.character(Biostrings::quality(Biostrings::quality(fq1))), as.character(Biostrings::quality(Biostrings::quality(fq2))))), 
      id = ShortRead::id(fq1))
    #Write to file
    dimsum__writeFastq(shortreads = fqc, outputFile = output_FASTQ, initial_write = initial_write)
    initial_write <- FALSE
    #Update statistics
    a_stats[['Merged']] <- a_stats[['Merged']] + length(fqc)
    a_stats[['merged_lengths']] <- c(a_stats[['merged_lengths']], IRanges::width(ShortRead::sread(fqc)))
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
}




