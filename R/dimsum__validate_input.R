
#' dimsum__validate_input
#'
#' Check validity and reformat input parameters.
#'
#' @param input_list a list of input arguments (required)
#'
#' @return an experiment metadata object
#' @export
dimsum__validate_input <- function(
  input_list
  ){
  #Check all input parameter types correct
  for(input_arg in names(input_list)){
    #Integer arguments
    if("integer" %in% input_list[[input_arg]][[2]]){
      #Check if numeric
      if(!typeof(input_list[[input_arg]][[1]]) %in% c(input_list[[input_arg]][[2]], "double")){
        stop(paste0("Invalid type of argument '", input_arg, "' (", input_list[[input_arg]][[2]], ")"), call. = FALSE)
      }else{
        #Check if whole number (if not NULL)
        if(!is.null(input_list[[input_arg]][[1]])){
          if(input_list[[input_arg]][[1]]%%1!=0){
            stop(paste0("Invalid type of argument '", input_arg, "' (", input_list[[input_arg]][[2]], ")"), call. = FALSE)
          }
          #Convert to integer
          input_list[[input_arg]][[1]] <- as.integer(input_list[[input_arg]][[1]])
        }
      }
    }else{
      #Other arguments
      if(!typeof(input_list[[input_arg]][[1]]) %in% input_list[[input_arg]][[2]]){
        stop(paste0("Invalid type of argument '", input_arg, "' (", input_list[[input_arg]][[2]], ")"), call. = FALSE)
      }
    }
  }
  #Metadata object
  dimsum_meta <- sapply(input_list, '[', 1)

  #Reformat directory paths and check if they exist (remove trailing "/" if present)
  for(input_arg in c("fastqFileDir", "outputPath")){
    dimsum_meta[[input_arg]] <- gsub("/$", "", dimsum_meta[[input_arg]])
    #Check if directory paths exist
    if(!file.exists(dimsum_meta[[input_arg]])){
      stop(paste0("Invalid '", input_arg, "' argument (directory not found)"), call. = FALSE)
    }
  }

  #FASTQ file extension
  #Append "." if not first character
  dimsum_meta[["fastqFileExtension"]] <- ifelse(
    substr(dimsum_meta[["fastqFileExtension"]], 1, 1)!=".", 
    paste0(".", dimsum_meta[["fastqFileExtension"]]), 
    dimsum_meta[["fastqFileExtension"]])
  #Check if alphanumeric
  if(!grepl("^[A-Za-z0-9]+$", substr(dimsum_meta[["fastqFileExtension"]], 2, nchar(dimsum_meta[["fastqFileExtension"]])), perl = T)){
    stop(paste0("Invalid '", "fastqFileExtension", "' argument (only alphanumeric characters allowed)"), call. = FALSE)
  }

  #Check WT sequence is valid (ACGT characters only) and save case-coded sequence
  all_characters <- unique(unlist(strsplit(dimsum_meta[["wildtypeSequence"]], "")))
  if(sum(!all_characters %in% c("A", "a", "C", "c", "G", "g", "T", "t"))!=0){
    stop("Invalid wild-type nucleotide sequence. Only valid nucleotide sequences allowed (A/C/G/T). Use lower-case letters (a/c/g/t) to indicate internal constant regions (removed before fitness calculations).", call. = FALSE)
  }else{
    #case-coded WT nucleotide sequence
    dimsum_meta[["wildtypeSequenceCoded"]] <- dimsum_meta[["wildtypeSequence"]]
    dimsum_meta[["wildtypeSequence"]] <- toupper(dimsum_meta[["wildtypeSequence"]])
  }

  #Check strictly positive integer usearch... arguments
  if(sum(unlist(dimsum_meta[c("usearchMinQual", "usearchMinlen", "usearchMinovlen")])<=0)!=0){
    stop("Invalid 'usearch...' arguments. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check positive integer fitness... arguments
  if(sum(unlist(dimsum_meta[c("fitnessMinInputCountAll", "fitnessMinInputCountAny", "fitnessHighConfidenceCount", "fitnessDoubleHighConfidenceCount")])<0)!=0){
    stop("Invalid 'fitness...Count...' arguments. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  #Check strictly positive integer splitChunkSize argument
  if(dimsum_meta[["splitChunkSize"]]<=0){
    stop("Invalid 'splitChunkSize' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check fitnessMaxSubstitutions argument greater than 1
  if(dimsum_meta[["fitnessMaxSubstitutions"]]<2){
    stop("Invalid 'fitnessMaxSubstitutions' argument. Only integers greater than 1 allowed.", call. = FALSE)
  }

  #Check errorModelMaxSubstitutions argument greater than 0
  if(dimsum_meta[["errorModelMaxSubstitutions"]]<1){
    stop("Invalid 'errorModelMaxSubstitutions' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check usearchMaxee argument
  if(dimsum_meta[["usearchMaxee"]]<=0){
    stop("Invalid 'usearchMaxee' argument. Only positive doubles allowed (zero exclusive).", call. = FALSE)
  }

  #Disable bayesianDoubleFitness option (TEMPORARY FIX: still in development)
  dimsum_meta[["bayesianDoubleFitness"]] <- FALSE

  #Check bayesianDoubleFitnessLamD argument
  if(dimsum_meta[["bayesianDoubleFitnessLamD"]]<=0){
    stop("Invalid 'bayesianDoubleFitnessLamD' argument. Only positive doubles allowed (zero exclusive).", call. = FALSE)
  }

  #Check startStage argument
  if(dimsum_meta[["startStage"]]<=0){
    stop("Invalid 'startStage' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check stopStage argument
  if(dimsum_meta[["stopStage"]]<0){
    stop("Invalid 'stopStage' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  #Check numCores argument
  if(dimsum_meta[["numCores"]]<=0){
    stop("Invalid 'numCores' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check library design paired if trans library specified
  if(dimsum_meta[["transLibrary"]] & !dimsum_meta[["paired"]]){
    stop("Invalid 'paired' argument. Only paired-end trans libraries allowed.", call. = FALSE)
  }

  #Check library design paired if unstranded library specified
  if(!dimsum_meta[["stranded"]] & !dimsum_meta[["paired"]]){
    stop("Invalid 'stranded' argument. Only unstranded paired-end libraries allowed.", call. = FALSE)
  }

  #Check sequenceType one of noncoding/coding/auto
  if(!dimsum_meta[["sequenceType"]] %in% c("nonconding", "coding", "auto")){
    stop("Invalid 'sequenceType' argument. Only noncoding/coding/auto allowed.", call. = FALSE)
  }else if(dimsum_meta[["sequenceType"]]=="auto"){
    dimsum_meta[["sequenceType"]] <- dimsum__detect_sequence_type(dimsum_meta[["wildtypeSequence"]])
  }

  #Check retainedReplicates comma-separated list of integers or "all"
  if(dimsum_meta[["retainedReplicates"]] != "all"){
    if(sum(!strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]] %in% strsplit("0123456789", "")[[1]])!=0){
      stop("Invalid 'retainedReplicates' argument. Only comma-separated list of integers or 'all' allowed.", call. = FALSE)
    }
  }

  #Return
  return(dimsum_meta)
}
