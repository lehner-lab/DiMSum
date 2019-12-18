
#' dimsum__check_barcode_design
#'
#' Validate metadata from barcode design file.
#'
#' @param barcode_design barcode design data.frame (required)
#'
#' @return Nothing
#' @export
dimsum__check_barcode_design <- function(
  barcode_design
  ){

  ### Column existence checks 
  #Check if mandatory columns present
  mandatory_cols <- c("pair1", "pair2", "barcode1", "barcode2", "new_pair_prefix")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(barcode_design)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from barcodeDesign file ('pair1', 'pair2', 'barcode1', 'barcode2', 'new_pair_prefix')"), call. = FALSE)
  }

  ### FASTQ file checks (pair1/pair2 columns)
  #Check file name columns are of type character 
  if(typeof(barcode_design[,"pair1"])!="character" | typeof(barcode_design[,"pair2"])!="character"){
    stop(paste0("One or more invalid FASTQ file name values in barcodeDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for incorrect pairing of FASTQ files
  unique_pairs <- unique(barcode_design[,c("pair1", "pair2")])
  if(sum(duplicated(unique_pairs[,"pair1"]))!=0 | sum(duplicated(unique_pairs[,"pair2"]))!=0){
    stop(paste0("FASTQ files not correctly paired in barcodeDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }

  ### Barcode checks
  #Check barcode column is of type character 
  if(typeof(barcode_design[,"barcode1"])!="character" | typeof(barcode_design[,"barcode2"])!="character"){
    stop(paste0("One or more invalid barcode values in barcodeDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check barcode sequences are valid (ACGT characters only)
  all_characters <- unique(unlist(strsplit(as.character(unlist(barcode_design[,c("barcode1", "barcode2")])), "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T", NA))!=0){
    stop("Invalid barcode sequences. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  ### New pair prefix name checks (new_pair_prefix column)
  #Check new_pair_prefix column is of type character 
  if(typeof(barcode_design[,"new_pair_prefix"])!="character"){
    stop(paste0("One or more invalid 'new_pair_prefix' values in barcodeDesign file (only alphanumeric and underscore characters allowed)"), call. = FALSE)
  }
  #Check for duplicated prefices
  if(sum(duplicated(barcode_design[,"new_pair_prefix"]))!=0){
    stop(paste0("Duplicate 'new_pair_prefix' values not allowed in barcodeDesign file"), call. = FALSE)
  }
  #Check only alphanumeric characters (and underscore)
  if(sum(!grepl("^[A-Za-z0-9_]+$", barcode_design[,"new_pair_prefix"], perl = T))!=0){
    stop(paste0("One or more invalid 'new_pair_prefix' values in barcodeDesign file (only alphanumeric and underscore characters allowed)"), call. = FALSE)
  }

  ### Duplicate matrix row checks
  #Check for duplicated rows in the following column pairs: "pair1":"barcode" and "pair2":"barcode"
  if(sum(duplicated(barcode_design[,c("pair1", "barcode1", "barcode2")]))!=0 | sum(duplicated(barcode_design[,c("pair2", "barcode1", "barcode2")]))!=0){
    stop(paste0("One or more duplicated mappings in barcodeDesign file matrix (FASTQ file to barcode mapping should be unique)"), call. = FALSE)
  }

}


