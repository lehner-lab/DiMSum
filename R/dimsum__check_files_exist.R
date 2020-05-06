
#' dimsum__check_files_exist
#'
#' Check files exist and if not produce error message and stop.
#'
#' @param required_files character vector of file paths (required)
#' @param stage_number integer stage number (required)
#' @param execute whether this stage is to be executed (default:TRUE)
#' @param exit whether to stop execution (default:TRUE)
#'
#' @return Nothing
#' @export
dimsum__check_files_exist <- function(
  required_files,
  stage_number = NULL,
  execute = TRUE,
  exit = TRUE
  ){
	missing_files <- required_files[!file.exists(required_files)]
  if(length(missing_files)!=0 & execute){
    if(is.null(stage_number)){
      #No stage number supplied
      dimsum__status_message(paste0("The following files do not exist.\n"))
      dimsum__status_message(paste0(missing_files, "\n"))
      if(exit){
        stop(paste0("Cannot proceed with DiMSum. Required files do not exist."), call. = FALSE)
      }
    }else{
      #Stage number supplied
      dimsum__status_message(paste0("The following files required for DiMSum STAGE ", stage_number, " do not exist.\n"))
      dimsum__status_message(paste0(missing_files, "\n"))
      dimsum__status_message(paste0("Ensure that the 'startStage' argument has been correctly specified and that all previous stages have successfully completed.\n"))
      dimsum__status_message(paste0("If you previously ran DiMSum with 'retainIntermediateFiles' set to FALSE, you will need set 'startStage' to 0 in order to regenerate the files required for this stage.\n"))
      if(exit){
        stop(paste0("Cannot proceed with DiMSum STAGE ", stage_number, ". Required files do not exist."), call. = FALSE)
      }
    }
  }
}
