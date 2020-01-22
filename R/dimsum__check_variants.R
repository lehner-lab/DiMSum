
#' dimsum__check_variants
#'
#' Check whether minimum variant requirements satisfied.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__check_variants <- function(
  dimsum_meta,
  input_dt
  ){

  #Check if WT sequence and at least 2 substitutions variants exist
  if(input_dt[WT==T,.N]==0){
    stop(paste0("Cannot proceed with fitness estimation: WT variant not found"), call. = FALSE)
  }else if(input_dt[is.na(WT),.N]<2){
    stop(paste0("Cannot proceed with fitness estimation: insufficient substitution variants"), call. = FALSE)
  }

}
