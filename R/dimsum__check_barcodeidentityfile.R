
#' dimsum__check_barcodeidentityfile
#'
#' Check whether user-specified barcode identity file correctly formatted.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param return_data whether or not to return data (default: FALSE)
#'
#' @return Reformatted data.table
#' @export
#' @import data.table
dimsum__check_barcodeidentityfile <- function(
  dimsum_meta,
  return_data = FALSE
  ){

  ### Abort if no count file supplied
  if(is.null(dimsum_meta[["barcodeIdentityPath"]])){return(NULL)}

  #Load file
  barcode_dt <- fread(dimsum_meta[["barcodeIdentityPath"]])
  #Check if mandatory columns present
  mandatory_cols <- c("barcode", "variant")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(barcode_dt)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from barcodeIdentity file ('barcode', 'variant')"), call. = FALSE)
  }
  #Check all columns are of type character 
  if(typeof(barcode_dt[,barcode])!="character" | typeof(barcode_dt[,variant])!="character"){
    stop(paste0("One or more invalid values in barcodeIdentity file (only A/C/G/T characters allowed)"), call. = FALSE)
  }

  ### Return data if required
  if(return_data){
    return(barcode_dt)
  }
}
