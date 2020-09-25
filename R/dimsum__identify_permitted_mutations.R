
#' dimsum__identify_permitted_mutations
#'
#' Identify variants with permitted mutations (if same length as WT).
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#'
#' @return A data.table with permitted variants identified
#' @export
#' @import data.table
dimsum__identify_permitted_mutations <- function(
  dimsum_meta,
  input_dt
  ){

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
  input_dt[indel==F, permitted := T]
  for(i in 1:nchar(dimsum_meta[["permittedSequences"]])){
    input_dt[indel==F & !toupper(substr(nt_seq, i, i)) %in% nuc_codes[[substr(dimsum_meta[["permittedSequences"]], i, i)]], permitted := F]
  }

  return(input_dt)

}
