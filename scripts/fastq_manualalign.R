
#fastq_manualalign
#
# Attempt to merge paired end reads according to a specified alignment length.
#
# input_FASTQ1: Path to first read FASTQ file (required)
# input_FASTQ2: Path to second read FASTQ file (required)
# output_FASTQ: Path to output FASTQ file (required)
# output_REPORT: Path to report file (required)
# num_nuc: Enforced alignment length in base pairs (required)
# min_qual: Minimum observed base quality to retain read pair (required)
# max_ee: Maximum number of expected errors to retain read pair (required)
# min_len: Discard pair if either read is shorter than this (required)
# concatentate_reads: Simply concatentate reads (without reverse complementing second read in pair) (default: FALSE)
#
# Returns: nothing.
#
fastq_manualalign <- function(
  input_FASTQ1,
  input_FASTQ2,
  output_FASTQ,
  output_REPORT,
  num_nuc,
  min_qual,
  max_ee,
  min_len,
  concatentate_reads = FALSE
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
  fq1 <- readFastq(input_FASTQ1)
  fq2 <- readFastq(input_FASTQ2)
  #Update statistics
  a_stats[['Pairs']] <- length(fq1)
  #Read lengths
  fq1_lengths <- width(sread(fq1))
  fq2_lengths <- width(sread(fq2))
  #Update statistics
  a_stats[['Fwd_too_short']] <- sum(fq1_lengths<min_len)
  a_stats[['Rev_too_short']] <- sum(fq2_lengths<min_len)
  #Subset to sequences at least 64 bp long
  fq1 <- fq1[fq1_lengths>=min_len & fq2_lengths>=min_len]
  fq2 <- fq2[fq1_lengths>=min_len & fq2_lengths>=min_len]
  #Simply paste reads
  if(concatentate_reads){
    #Update statistics
    a_stats[['Alignments_zero_diffs']] <- sum(fq1_lengths>=min_len & fq2_lengths>=min_len)
    #Read quality matrices
    qmat1 <- as(quality(fq1), "matrix")
    qmat2 <- as(quality(fq2), "matrix")
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
    fqc <- ShortReadQ(
      sread = DNAStringSet(paste0(as.character(sread(fq1)), as.character(sread(fq2)))), 
      quality = BStringSet(paste0(as.character(quality(quality(fq1))), as.character(quality(quality(fq2))))), 
      id = ShortRead::id(fq1))
    #Write to file
    writeFastq(fqc, file = output_FASTQ, compress = F)
    #Update statistics
    a_stats[['Merged']] = length(fqc)
    if(a_stats[['Merged']]!=0){
      a_stats[['Merged_length_min']] = min(width(sread(fqc)))
      a_stats[['Merged_length_low']] = as.numeric(quantile(width(sread(fqc)), 0.25))
      a_stats[['Merged_length_median']] = as.numeric(median(width(sread(fqc))))
      a_stats[['Merged_length_high']] = as.numeric(quantile(width(sread(fqc)), 0.75))
      a_stats[['Merged_length_max']] = max(width(sread(fqc)))
    }else{
      a_stats[['Merged_length_min']] = NA
      a_stats[['Merged_length_low']] = NA
      a_stats[['Merged_length_median']] = NA
      a_stats[['Merged_length_high']] = NA
      a_stats[['Merged_length_max']] = NA
    }
  }else{
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




