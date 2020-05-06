
#' dimsum__usearch_single_end_library_helper
#'
#' Helper function to filter single end reads. "dimsum_meta", "usearch_outpath", "sample_names" and "dimsum__filter_single_end_reads" objects need to be available globally.
#'
#' @param i "sample_names" index (required)
#'
#' @return Nothing
#' @export
dimsum__usearch_single_end_library_helper <- function(
  i
  ){
  temp_out <- dimsum__filter_single_end_reads(
    input_FASTQ = file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair1"]),
    output_FASTQ = file.path(usearch_outpath, paste0(sample_names[i], '.usearch')),
    output_REPORT = file.path(usearch_outpath, paste0(sample_names[i], '.report')),
    min_qual = dimsum_meta[["usearchMinQual"]],
    max_ee = dimsum_meta[["usearchMaxee"]],
    min_len = dimsum_meta[["usearchMinlen"]])
}
