
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

  #Check whether experiment design is valid
  dimsum__check_barcode_design(barcode_design)

  #Check FASTQ files exist
  #Pair1 files
  for(i in unlist(barcode_design[,c("pair1")])){
    if(!file.exists(file.path(dimsum_meta[["fastqFileDir"]], i))){
      stop(paste0("Invalid FASTQ file name '", i, "' in barcodeDesign file (file not found)"), call. = FALSE)
    }
  }
  #Pair2 files
  for(i in unlist(barcode_design[,c("pair2")])){
    if(!file.exists(file.path(dimsum_meta[["fastqFileDir"]], i))){
      stop(paste0("Invalid FASTQ file name '", i, "' in barcodeDesign file (file not found)"), call. = FALSE)
    }
  }

  return(barcode_design)
}


