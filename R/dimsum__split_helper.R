
#' dimsum__split_helper
#'
#' Helper function to split fastq files. "dimsum_meta", "fastq_pair_list", "split_outpath", "dimsum__fastq_splitter" and "dimsum__writeFastq" objects need to be available globally.
#'
#' @param i "fastq_pair_list" row index (required)
#'
#' @return Nothing
#' @export
dimsum__split_helper <- function(
  i
  ){
  num_records <- dimsum__fastq_splitter(
    inputFile = file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[i,][1]),
    outputFilePrefix = file.path(split_outpath, paste0(fastq_pair_list[i,][1], ".split")),
    chunkSize = dimsum_meta[["splitChunkSize"]])
  if(dimsum_meta[["paired"]]){
    num_records <- dimsum__fastq_splitter(
      inputFile = file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[i,][2]),
      outputFilePrefix = file.path(split_outpath, paste0(fastq_pair_list[i,][2], ".split")),
      numRecords = num_records)
  }
}
