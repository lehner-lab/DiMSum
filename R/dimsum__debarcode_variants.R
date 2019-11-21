
#' dimsum__debarcode_variants
#'
#' Debarcode variants.
#'
#' @param input_dt input data.table (required)
#' @param barcode_path path to barcode identity file (required)
#'
#' @return A data.table with variants debarcoded
#' @export
#' @import data.table
dimsum__debarcode_variants <- function(
  input_dt,
  barcode_path
  ){

  #Load barcode to variant mapping
  barcode_dt <- data.table::fread(barcode_path)
  barcode_list <- as.list(tolower(unlist(barcode_dt[,variant])))
  names(barcode_list) <- tolower(unlist(barcode_dt[,barcode]))
  
  #Variants with valid barcodes
  input_dt[!nt_seq %in% names(barcode_list), barcode_valid := FALSE]

  #Debarcode variants
  debarcode_dt <- input_dt[barcode_valid==T,]
  debarcode_dt[, nt_seq := unlist(barcode_list[nt_seq])]
  #Aggregate counts for identical variants
  idx <- names(debarcode_dt)[grep(names(debarcode_dt),pattern="_count$")]
  for (i in seq_along(idx)) {
    #Aggregate counts accross identical AA variants
    debarcode_dt[,paste0(idx[i],"_agg") := sum(.SD),nt_seq,.SDcols = idx[i]]
  }
  #Retain only one row per AA variant
  debarcode_dt <- debarcode_dt[!duplicated(nt_seq),.SD,nt_seq,.SDcols = c(
    names(debarcode_dt)[grep(names(debarcode_dt),pattern="_agg$")],
    "barcode_valid")]
  #Rename count columns
  names(debarcode_dt) <- gsub("_agg$", "", names(debarcode_dt))

  #Merge variant data
  output_dt <- rbind(debarcode_dt, input_dt[barcode_valid==F,])

  return(output_dt)

}
