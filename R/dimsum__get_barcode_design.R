
#' dimsum__get_barcode_design
#'
#' Get, format and validate metadata from barcode design file.
#'
#' @param dimsum_meta an experiment metadata object (required)
#'
#' @return a data.frame with the validated barcode design 
#' @export
dimsum__get_barcode_design <- function(
  dimsum_meta
  ){
  #Check if barcodeDesignPath specified
  if(is.null(dimsum_meta[["barcodeDesignPath"]])){
    return(NULL)
  }

  #Load barcode design
  if(!file.exists(dimsum_meta[["barcodeDesignPath"]])){
    stop(paste0("Invalid '", "barcodeDesignPath", "' argument (file not found)"), call. = FALSE)
  }
  barcode_design <- read.table(dimsum_meta[["barcodeDesignPath"]], header = T, stringsAsFactors = F, sep="\t")

  #Set pair2 column equal to pair1 column if single-end library (contents of existing pair2 column will be ignored)
  if(!dimsum_meta[["paired"]]){
    if(!"pair1" %in% colnames(barcode_design)){
      stop(paste0("Mandatory column missing from barcodeDesign file ('pair1')"), call. = FALSE)
    }else{
      barcode_design[,"pair2"] <- barcode_design[,"pair1"]
    }
  }

  #Add original FASTQ directory
  if("fastqFileDir" %in% names(barcode_design)){
    names(barcode_design)[names(barcode_design)=="fastqFileDir"] <- "pair_directory"
  }else{
    barcode_design[,"pair_directory"] <- dimsum_meta[["fastqFileDir"]]
  }

  #Set barcode1 column equal to barcode column (backwards compatibility)
  if("barcode" %in% colnames(barcode_design)){
    barcode_design[,"barcode1"] <- barcode_design[,"barcode"]
  }

  #Set barcode1 column equal to barcode column (backwards compatibility)
  if(!"barcode2" %in% colnames(barcode_design)){
    barcode_design[,"barcode2"] <- barcode_design[,"barcode1"]
  }

  #Check whether barcode design is valid
  dimsum__check_barcode_design(barcode_design)

  #Check FASTQ files exist
  #Pair1 files
  for(i in unlist(file.path(barcode_design[,"pair_directory"], barcode_design[,"pair1"]))){
    if(!file.exists(i)){
      stop(paste0("Invalid FASTQ file name '", i, "' in barcodeDesign file (file not found)"), call. = FALSE)
    }
  }
  #Pair2 files
  for(i in unlist(file.path(barcode_design[,"pair_directory"], barcode_design[,"pair2"]))){
    if(!file.exists(i)){
      stop(paste0("Invalid FASTQ file name '", i, "' in barcodeDesign file (file not found)"), call. = FALSE)
    }
  }

  return(barcode_design)
}


