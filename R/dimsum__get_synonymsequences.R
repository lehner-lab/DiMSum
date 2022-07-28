
#' dimsum__get_synonymsequences
#'
#' Get and check whether user-specified synonym sequences correctly formatted.
#'
#' @param dimsum_meta an experiment metadata object (required)
#'
#' @return character vector of sequences
#' @export
#' @import data.table
dimsum__get_synonymsequences <- function(
  dimsum_meta
  ){

  ### Abort if no count file supplied
  if(is.null(dimsum_meta[["synonymSequencePath"]])){return(NULL)}

  input_dt <- data.table::fread(dimsum_meta[["synonymSequencePath"]], header = F)

  ### Nucleotide sequence checks (single column)
  #Check if just a single columns present
  if(nrow(input_dt)==0 | ncol(input_dt)!=1){
    stop(paste0("Incorrect number of columns in file specified by synonymSequencePath (only one allowed)"), call. = FALSE)
  }
  #Check nucleotide sequence column is of type character 
  if(sapply(input_dt, typeof)["V1"]!="character"){
    stop("One or more invalid sequences in synonym sequences file specified by 'synonymSequencePath'. Only valid coding nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }
  #Set nucleotide sequence to lower case
  input_dt[, V1 := tolower(V1)]
  #Check nucleotide sequences are valid (ACGT characters only)
  if(sum(!input_dt[,unique(unlist(strsplit(V1, "")))] %in% c('a', 'c', 'g', 't'))!=0){
    stop("One or more invalid sequences in synonym sequences file specified by 'synonymSequencePath'. Only valid coding nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  #Check if sequences are coding
  if(sum(unlist(lapply(as.list(input_dt[,V1]), dimsum__detect_sequence_type))!='coding')!=0){
    stop("One or more non-coding sequences in synonym sequences file specified by 'synonymSequencePath'. Only valid coding nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  #Translate sequences
  suppressWarnings(input_dt[, aa_seq := as.character(Biostrings::translate(Biostrings::DNAStringSet(V1), no.init.codon=T))])

  return(input_dt[,aa_seq])
}
