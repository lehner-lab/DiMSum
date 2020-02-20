
#' dimsum__convert_linked_adapters
#'
#' Convert to linked adapters if 3' constant region expected to be sequenced.
#'
#' @param dimsum_meta an experiment metadata object (required)
#'
#' @return An updated dimsum_meta object
#' @export
dimsum__convert_linked_adapters <- function(
  dimsum_meta){
  #Loop over all experimental design entries
  for(i in 1:dim(dimsum_meta[['exp_design']])[1]){
    num_cut5f <- ifelse(is.na(dimsum_meta[['exp_design']][i,"cutadaptCut5First"]), 0, dimsum_meta[['exp_design']][i,"cutadaptCut5First"])
    num_cut5s <- ifelse(is.na(dimsum_meta[['exp_design']][i,"cutadaptCut5Second"]), 0, dimsum_meta[['exp_design']][i,"cutadaptCut5Second"])
    num_cut3f <- ifelse(is.na(dimsum_meta[['exp_design']][i,"cutadaptCut3First"]), 0, dimsum_meta[['exp_design']][i,"cutadaptCut3First"])
    num_cut3s <- ifelse(is.na(dimsum_meta[['exp_design']][i,"cutadaptCut3Second"]), 0, dimsum_meta[['exp_design']][i,"cutadaptCut3Second"])
    #Read 1
    #Check if 5' constant region specified
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt5First"]) & !grepl("\\.\\.\\.", dimsum_meta[['exp_design']][i,"cutadapt5First"]) ){
      if( (dimsum_meta[['exp_design']][i,"pair1_length"]-num_cut5f-num_cut3f) > (nchar(dimsum_meta[['exp_design']][i,"cutadapt5First"]) + nchar(dimsum_meta[['wildtypeSequence']])) ){
        dimsum_meta[['exp_design']][i,"cutadapt5First"] <- paste0(dimsum_meta[['exp_design']][i,"cutadapt5First"], "...", dimsum_meta[['exp_design']][i,"cutadapt3First"])
        dimsum_meta[['exp_design']][i,"cutadapt3First"] <- NA
      }
    }
    #Read2 (read1 lengths can be variable due to inconsistent barcode trimming with cutadapt version <2.3)
    #Check if 5' constant region specified
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt5Second"]) & !grepl("\\.\\.\\.", dimsum_meta[['exp_design']][i,"cutadapt5Second"]) ){
      if( (dimsum_meta[['exp_design']][i,"pair2_length"]-num_cut5s-num_cut3s) > (nchar(dimsum_meta[['exp_design']][i,"cutadapt5Second"]) + nchar(dimsum_meta[['wildtypeSequence']])) ){
        dimsum_meta[['exp_design']][i,"cutadapt5Second"] <- paste0(dimsum_meta[['exp_design']][i,"cutadapt5Second"], "...", dimsum_meta[['exp_design']][i,"cutadapt3Second"])
        dimsum_meta[['exp_design']][i,"cutadapt3Second"] <- NA
      }
    }
  }
  return(dimsum_meta)
}
