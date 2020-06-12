
#' dimsum__unique_helper
#'
#' Helper function to concatenate split files. "dimsum_meta" and "unique_outpath" objects need to be available globally.
#'
#' @param i "sample_names" index (required)
#'
#' @return Nothing
#' @export
dimsum__unique_helper <- function(
  i
  ){
  this_sample_code <- unique(dimsum_meta[["exp_design"]][,"sample_code"])[i]
  read_pairs <- dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"sample_code"]==this_sample_code,"aligned_pair"]
  #Concatenate FASTQ files
  output_file1 <- gsub("_split1.vsearch$", ".vsearch", read_pairs[1])
  file.copy(
    from = file.path(dimsum_meta[["exp_design"]][,"aligned_pair_directory"][1], read_pairs[1]),
    to = file.path(unique_outpath, output_file1))
  if(length(read_pairs)>1){
    for(j in 2:length(read_pairs)){
      file.append(
        file1 = file.path(unique_outpath, output_file1),
        file2 = file.path(dimsum_meta[["exp_design"]][,"aligned_pair_directory"][1], read_pairs[j]))
    }
  }
}
