
#' dimsum__vsearch_trans_library_helper
#'
#' Helper function to merge trans library reads. "dimsum_meta", "vsearch_outpath", "sample_names" and "dimsum__concatenate_reads" objects need to be available globally.
#'
#' @param i "sample_names" index (required)
#'
#' @return Nothing
#' @export
dimsum__vsearch_trans_library_helper <- function(
  i
  ){
  temp_out <- dimsum__concatenate_reads(
    input_FASTQ1 = file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair1"]),
    input_FASTQ2 = file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair2"]),
    output_FASTQ = file.path(vsearch_outpath, paste0(sample_names[i], '.vsearch.gz')),
    output_REPORT = file.path(vsearch_outpath, paste0(sample_names[i], '.report')),
    min_qual = dimsum_meta[["vsearchMinQual"]],
    max_ee = dimsum_meta[["vsearchMaxee"]],
    min_len = dimsum_meta[["cutadaptMinLength"]],
    reverse_complement_second_read = dimsum_meta[["transLibraryReverseComplement"]])
}
