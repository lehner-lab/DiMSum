
#' save_metadata
#'
#' Save experiment workspace including metadata and session information.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param n the number of generations to go back (required)
#'
#' @return Nothing
#' @export
save_metadata <- function(
  dimsum_meta,
  n
  ){
  save(
    list = c(ls(parent.frame(n = n))), 
    file = file.path(dimsum_meta[["project_path"]], paste0(dimsum_meta[["project_name"]], '_workspace.RData')),
    envir = parent.frame(n = n))
  session_info = sessionInfo()
  save(session_info, file = file.path(dimsum_meta[["project_path"]], paste0(dimsum_meta[["project_name"]], '_sessionInfo.RData')))
}
