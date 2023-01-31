
#' dimsum__filter_reads_helper
#'
#' Helper function to filter reads. "dimsum_meta", "vsearch_outpath", "sample_names" and "dimsum__filter_reads" objects need to be available globally.
#'
#' @param i "sample_names" index (required)
#'
#' @return Nothing
#' @export
dimsum__filter_reads_helper <- function(
  i
  ){
  temp_out <- dimsum__filter_reads(
    input_FASTQ = file.path(vsearch_outpath, paste0(sample_names[i], '.vsearch.prefilter.gz')),
    input_REPORT = file.path(vsearch_outpath, paste0(sample_names[i], '.report.prefilter')),
    output_FASTQ = file.path(vsearch_outpath, paste0(sample_names[i], '.vsearch.gz')),
    output_REPORT = file.path(vsearch_outpath, paste0(sample_names[i], '.report')),
    min_qual = dimsum_meta[["vsearchMinQual"]])
}
