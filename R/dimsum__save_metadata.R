
#' dimsum__save_metadata
#'
#' Save experiment workspace including metadata and session information.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param n the number of generations to go back (required)
#'
#' @return Nothing
#' @export
dimsum__save_metadata <- function(
  dimsum_meta,
  n
  ){
  save(
    list = c(ls(parent.frame(n = n))), 
    file = file.path(dimsum_meta[["project_path"]], paste0(dimsum_meta[["projectName"]], '_workspace.RData')),
    envir = parent.frame(n = n),
    version = 2)
  session_info <- sessionInfo()
  save(session_info, 
  	file = file.path(dimsum_meta[["project_path"]], paste0(dimsum_meta[["projectName"]], '_sessionInfo.RData')),
  	version = 2)
}
