
#' dimsum__plot_samplename
#'
#' Format sample names for plotting purposes.
#'
#' @param input_names a list of sample names (required)
#'
#' @return reformatted sample names
#' @export
dimsum__plot_samplename <- function(
  input_names){
  if(sum(grepl("_t", input_names)) != 0){
    output_names <- sapply(lapply(lapply(strsplit(input_names, "_"), unlist), '[', c(1,5)), paste, collapse = "_")
    output_names <- gsub("_tNA", "", output_names)
  }else{
    output_names <- sapply(lapply(strsplit(input_names, "_"), unlist), '[', 1)
  }
  return(output_names)
}
