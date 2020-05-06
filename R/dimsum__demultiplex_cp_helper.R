
#' dimsum__demultiplex_cp_helper
#'
#' Helper function to copy and rename fastq files before demultiplexing. "dimsum_meta","fastq_pair_list" and "demultiplex_outpath" objects need to be available globally.
#'
#' @param i "fastq_pair_list" row index (required)
#'
#' @return Nothing
#' @export
dimsum__demultiplex_cp_helper <- function(
  i
  ){
  pair_name <- rownames(fastq_pair_list)[i]
  temp_design <- dimsum_meta[['barcode_design']][dimsum_meta[['barcode_design']][,'pair1']==fastq_pair_list[pair_name,'pair1'],]
  write(
    x = c(rbind(paste0('>', temp_design[,"new_pair_prefix"]), paste0('^', temp_design[,"barcode1"]))), 
    file = file.path(demultiplex_outpath, paste0('demultiplex_barcode1-file_', pair_name, '.fasta')), 
    sep="\n")
  write(
    x = c(rbind(paste0('>', temp_design[,"new_pair_prefix"]), paste0('^', temp_design[,"barcode2"]))), 
    file = file.path(demultiplex_outpath, paste0('demultiplex_barcode2-file_', pair_name, '.fasta')), 
    sep="\n")
  #Check if file extension incompatible with cutadapt (i.e. NOT ".fastq" or ".fastq.gz")
  if(dimsum_meta[["fastqFileExtension"]]!=".fastq"){
    #Copy FASTQ files to temp directory and format extension
    new_fastq_name1 <- gsub(paste0(dimsum_meta[["fastqFileExtension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), c(".fastq", ".fastq.gz")[as.numeric(dimsum_meta[["gzipped"]])+1], fastq_pair_list[pair_name,][1])
    new_fastq_name2 <- gsub(paste0(dimsum_meta[["fastqFileExtension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), c(".fastq", ".fastq.gz")[as.numeric(dimsum_meta[["gzipped"]])+1], fastq_pair_list[pair_name,][2])
    file.copy(
      from = file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][1]),
      to = file.path(demultiplex_outpath, new_fastq_name1))
    #If second read in pair exists
    if(dimsum_meta[["paired"]]){
      file.copy(
        from = file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][2]),
        to = file.path(demultiplex_outpath, new_fastq_name2))
    }
  }
}
