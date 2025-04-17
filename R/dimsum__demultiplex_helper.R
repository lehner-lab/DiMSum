
#' dimsum__demultiplex_helper
#'
#' Helper function to demultiplex fastq files. "dimsum_meta", "demultiplex_outpath" and "fastq_pair_list" objects need to be available globally.
#'
#' @param i "fastq_pair_list" row index (required)
#'
#' @return Nothing
#' @export
dimsum__demultiplex_helper <- function(
  i
  ){
  pair_name <- rownames(fastq_pair_list)[i]
  #Demultiplex using cutadapt
  if(dimsum_meta[["paired"]]){
    temp_out <- system(paste0(
      "cutadapt",
      " -g file:",
      file.path(demultiplex_outpath, paste0('demultiplex_barcode1-file_', pair_name, '.fasta')),
      " -G file:",
      file.path(demultiplex_outpath, paste0('demultiplex_barcode2-file_', pair_name, '.fasta')),
      " -e ",
      as.character(dimsum_meta[["barcodeErrorRate"]]),
      " --no-indels ",
      " --pair-adapters ",
      " --untrimmed-output ",
      file.path(demultiplex_outpath, paste0(basename(fastq_pair_list[pair_name,][1]), ".demultiplex.unknown.fastq.gz")),
      " --untrimmed-paired-output ",
      file.path(demultiplex_outpath, paste0(basename(fastq_pair_list[pair_name,][2]), ".demultiplex.unknown.fastq.gz")),
      " -o ",
      file.path(demultiplex_outpath, "{name}1.fastq.gz"),
      " -p ",
      file.path(demultiplex_outpath, "{name}2.fastq.gz"),
      " ",
      fastq_pair_list[pair_name,][1],
      " ",
      fastq_pair_list[pair_name,][2],
      " > ",
      file.path(demultiplex_outpath, paste0(basename(fastq_pair_list[pair_name,][1]), ".demultiplex.stdout")),
      " 2> ",
      file.path(demultiplex_outpath, paste0(basename(fastq_pair_list[pair_name,][1]), ".demultiplex.stderr"))))
  }else{
    temp_out <- system(paste0(
      "cutadapt",
      " -g file:",
      file.path(demultiplex_outpath, paste0('demultiplex_barcode1-file_', pair_name, '.fasta')),
      " -e ",
      as.character(dimsum_meta[["barcodeErrorRate"]]),
      " --no-indels ",
      " --untrimmed-output ",
      file.path(demultiplex_outpath, paste0(basename(fastq_pair_list[pair_name,][1]), ".demultiplex.unknown.fastq.gz")),
      " -o ",
      file.path(demultiplex_outpath, "{name}1.fastq.gz"),
      " ",
      fastq_pair_list[pair_name,][1],
      " > ",
      file.path(demultiplex_outpath, paste0(basename(fastq_pair_list[pair_name,][1]), ".demultiplex.stdout")),
      " 2> ",
      file.path(demultiplex_outpath, paste0(basename(fastq_pair_list[pair_name,][1]), ".demultiplex.stderr"))))        
  }
}
