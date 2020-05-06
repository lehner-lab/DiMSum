
#' dimsum__unzip_helper
#'
#' Helper function to unzip fastq files. "all_fastq" and "fastq_outpath" objects need to be available globally.
#'
#' @param i "all_fastq" index (required)
#'
#' @return Nothing
#' @export
dimsum__unzip_helper <- function(
  i
  ){
  temp_out <- system(paste0(
    "gunzip -c ", 
    all_fastq[i], 
    " > ", 
    file.path(fastq_outpath, gsub(".gz$", "", basename(all_fastq[i]))),
    " 2> ",
    file.path(fastq_outpath, paste0(gsub(".gz$", "", basename(all_fastq[i])), '.stderr'))))
}
