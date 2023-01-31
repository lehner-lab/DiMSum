
#' dimsum__vsearch_single_end_library_helper
#'
#' Helper function to filter single end reads. "dimsum_meta", "vsearch_outpath", "sample_names" and "dimsum__filter_single_end_reads" objects need to be available globally.
#'
#' @param i "sample_names" index (required)
#'
#' @return Nothing
#' @export
dimsum__vsearch_single_end_library_helper <- function(
  i
  ){
  temp_out <- dimsum__filter_single_end_reads(
    input_FASTQ = file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair1"]),
    output_FASTQ = file.path(vsearch_outpath, paste0(sample_names[i], '.vsearch.gz')),
    output_REPORT = file.path(vsearch_outpath, paste0(sample_names[i], '.report')),
    min_qual = dimsum_meta[["vsearchMinQual"]],
    max_ee = dimsum_meta[["vsearchMaxee"]],
    min_len = dimsum_meta[["cutadaptMinLength"]])
}
