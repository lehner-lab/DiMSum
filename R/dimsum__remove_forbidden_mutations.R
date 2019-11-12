
#' dimsum__remove_forbidden_mutations
#'
#' Subset variants to those with permitted mutations.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#'
#' @return A data.table with variants with permitted mutations only
#' @export
#' @import data.table
dimsum__remove_forbidden_mutations <- function(
  dimsum_meta,
  input_dt
  ){

  message("Removing nucleotide variants without permitted mutations...")

  nuc_codes <- list(
    "A" = "A",
    "C" = "C",
    "G" = "G",
    "T" = "T",
    "R" = c("A", "G"),
    "Y" = c("C", "T"),
    "S" = c("C", "G"),
    "W" = c("A", "T"),
    "K" = c("G", "T"),
    "M" = c("A", "C"),
    "B" = c("C", "G", "T"),
    "D" = c("A", "G", "T"),
    "H" = c("A", "C", "T"),
    "V" = c("A", "C", "G"),
    "N" = c("A", "C", "G", "T")
    )

  #Loop over all nucleotide sequence positions and subset to permitted mutations
  for(i in 1:nchar(dimsum_meta[["permittedSequences"]])){
    input_dt <- input_dt[toupper(substr(nt_seq, i, i)) %in% nuc_codes[[substr(dimsum_meta[["permittedSequences"]], i, i)]]]
  }

  message("Done")

  #Output data.table
  output_dt <- input_dt

  return(output_dt)

}
