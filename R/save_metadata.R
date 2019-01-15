
#' save_metadata
#'
#' Save experiment metadata.
#'
#' @param dimsum_meta an experiment metadata object (required)
#'
#' @return Nothing
#' @export
save_metadata <- function(dimsum_meta){
  save(dimsum_meta, file=file.path(dimsum_meta[["project_path"]], paste0(dimsum_meta[["project_name"]], '_metadata.RData')))
}
