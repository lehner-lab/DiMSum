
#create_dimsum_dir
#
# Create results folder for dimsum pipeline stage.
#
# dimsum_dir: directory path string (required)
# execute: whether or not the system command will be executed (required)
# message: message string (optional, default: NULL i.e. no message displayed)
# overwrite_dir: delete directory if already exists (optional, default: TRUE)
#
# Returns: nothing.
#
create_dimsum_dir <- function(dimsum_dir, execute, message = NULL, overwrite_dir = TRUE){
  if(execute){
    if(!is.null(message)){
      message(paste("\n\n\n*******", message, "*******\n\n\n"))
    }
    if(dir.exists(dimsum_dir) & overwrite_dir){
      unlink(dimsum_dir, recursive = TRUE)
    }
    dir.create(dimsum_dir, showWarnings = FALSE)
  }else{
    if(!is.null(message)){
      message(paste("\n\n\n*******", message, "(not executed)", "*******\n\n\n"))
    }
  }
}
