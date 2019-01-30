
#' dimsum__fastq_splitter
#'
#' Split FASTQ file into roughly equally sized chunks (bytes).
#'
#' @param inputFile Path to FASTQ file (required)
#' @param outputFilePrefix Prefix to output FASTQ file (required)
#' @param chunkSize Chunk size in bytes
#' @param numRecords Number of FASTQ records per file
#'
#' @return Number of FASTQ records per file
#' @export
dimsum__fastq_splitter <- function(
  inputFile,
  outputFilePrefix,
  chunkSize=NULL,
  numRecords=NULL
  ){
  #Determine number of fastq records to write per output file
  num_fastqrecords <- NULL
  if(!is.null(numRecords)){
    #Number of records provided
    num_fastqrecords <- numRecords
  }else if(!is.null(chunkSize)){
    #Number of bytes per file provided
    #Open FASTQ file
    f <- ShortRead::FastqStreamer(inputFile, 1)
    fq <- ShortRead::yield(f)
    #Number of bytes in first FASTQ record
    bytesperrecord <- IRanges::width(ShortRead::id(fq)) + 1 + IRanges::width(ShortRead::sread(fq))*2 + 1 + 4
    #Estimated lines per file
    num_fastqrecords <- as.integer(chunkSize/bytesperrecord)
    close(f)
  }else{
    stop("Either --chunkSize or --numRecords arguments need to be specified", call. = FALSE)
  }

  #Split all FASTQ records into numbered output files with max num_fastqrecords per file
  count <- 1
  records <- 0 #records written to this file already
  output_FASTQ <- paste0(outputFilePrefix, count, '.fastq')
  #Open FASTQ file
  yield_size <- 1e6
  f <- ShortRead::FastqStreamer(inputFile, n=yield_size)
  while(length(fq <- ShortRead::yield(f))){
    #Write all yielded records
    if((records + length(fq)) >= num_fastqrecords){
      #File will be full after this write
      num_to_write <- num_fastqrecords-records
      #Write a subset of yielded records
      dimsum__writeFastq(shortreads = fq[1:num_to_write], outputFile = output_FASTQ, initial_write = (records==0))
      #Remove written records
      fq <- fq[-(1:num_to_write)]
      #Increment file counter and update output file
      count <- count + 1
      output_FASTQ <- paste0(outputFilePrefix, count, '.fastq')
      #Reset record counter
      records <- 0
    }else{
      #File will not be full after this write (write all yielded records)
      dimsum__writeFastq(shortreads = fq, outputFile = output_FASTQ, initial_write = (records==0))
      #Increment record counter
      records <- records + length(fq)
      #Remove written records
      fq <- fq[-(1:length(fq))]
    }
  }
  close(f)
  return(num_fastqrecords)
}

