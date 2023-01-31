
#' dimsum__writeFastq
#'
#' writeFastq wrapper.
#'
#' @param shortreads ShortReadQ object (required)
#' @param outputFile Path to output FASTQ file (required)
#' @param initial_write Is this the first write? (required)
#'
#' @return Nothing
#' @export
dimsum__writeFastq <- function(
  shortreads,
  outputFile,
  initial_write
  ){
  if(initial_write){
    if(file.exists(outputFile)){
      stop("file output path exists", call. = FALSE)
    }else{
      ShortRead::writeFastq(shortreads, outputFile, mode="w", compress=TRUE, full = FALSE)
    }
  }else{
    ShortRead::writeFastq(shortreads, outputFile, mode="a", compress=TRUE, full = FALSE)
  }
}
