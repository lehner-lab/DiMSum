
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

  input_list_original <- input_list

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

  #Save original command-line arguments
  dimsum_meta[["arg_list"]] <- sapply(input_list_original, '[', 1)

  #Either fastqFileDir or countPath need to be specified
  if(!is.null(dimsum_meta[["fastqFileDir"]])){
    #Reformat fastqFileDir and check if exists (remove trailing "/" if present)
    dimsum_meta[["fastqFileDir"]] <- gsub("/$", "", dimsum_meta[["fastqFileDir"]])
    if(!file.exists(dimsum_meta[["fastqFileDir"]])){
      stop(paste0("Invalid 'fastqFileDir' argument (directory not found)"), call. = FALSE)
    }
  }else if(!is.null(dimsum_meta[["countPath"]])){
    #Check if exists
    if(!file.exists(dimsum_meta[["countPath"]])){
      stop(paste0("Invalid 'countPath' argument (file not found)"), call. = FALSE)
    }
  }else{
    stop(paste0("Either 'fastqFileDir' or 'countPath' arguments must be specified"), call. = FALSE)
  }

  #Check if synonymSequencePath specified
  if(!is.null(dimsum_meta[["synonymSequencePath"]])){
    #Check if exists
    if(!file.exists(dimsum_meta[["synonymSequencePath"]])){
      stop(paste0("Invalid 'synonymSequencePath' argument (file not found)"), call. = FALSE)
    }
  }

  #Reformat outputPath and check if exists (remove trailing "/" if present)
  dimsum_meta[["outputPath"]] <- gsub("/$", "", dimsum_meta[["outputPath"]])
  if(!file.exists(dimsum_meta[["outputPath"]])){
    stop(paste0("Invalid 'outputPath' argument (directory not found)"), call. = FALSE)
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

  #WT sequence
  #Check WT sequence is valid (ACGT characters only) and save case-coded sequence
  all_characters <- unique(unlist(strsplit(dimsum_meta[["wildtypeSequence"]], "")))
  if(sum(!all_characters %in% c("A", "a", "C", "c", "G", "g", "T", "t"))!=0){
    stop("Invalid wild-type nucleotide sequence. Only valid nucleotide sequences allowed (A/C/G/T). Use lower-case letters (a/c/g/t) to indicate internal constant regions (removed before fitness calculations).", call. = FALSE)
  }
  #Check WT sequence contains at least one mutated position
  if(sum(all_characters %in% c("A", "C", "G", "T"))==0){
    stop("Invalid wild-type nucleotide sequence. Must contain at least one mutated position (A/C/G/T).", call. = FALSE)
  }
  #Case-coded WT nucleotide sequence
  dimsum_meta[["wildtypeSequenceCoded"]] <- dimsum_meta[["wildtypeSequence"]]
  dimsum_meta[["wildtypeSequence"]] <- toupper(dimsum_meta[["wildtypeSequence"]])

  #permittedSequences
  all_characters <- unlist(strsplit(dimsum_meta[["wildtypeSequenceCoded"]], ""))
  #Set permittedSequences argument if not specified
  if(is.null(dimsum_meta[["permittedSequences"]])){
    dimsum_meta[["permittedSequences"]] <- paste0(rep("N", sum(all_characters %in% c("A", "C", "G", "T"))), collapse = "")
  }
  #Check permittedSequences argument has same length as mutated sequence
  if(nchar(dimsum_meta[["permittedSequences"]]) != sum(all_characters %in% c("A", "C", "G", "T"))){
    stop("Invalid 'permittedSequences' argument length. Length should match the number of mutated positions in the wild-type nucleotide sequence (A/C/G/T).", call. = FALSE)
  }
  #Check permittedSequences argument sequence is valid
  all_characters <- unique(unlist(strsplit(dimsum_meta[["permittedSequences"]], "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"))!=0){
    stop("Invalid 'permittedSequences' argument. Only valid nucleotide codes allowed (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N).", call. = FALSE)
  }

  #Check strictly positive integer vsearch... arguments
  if(sum(unlist(dimsum_meta[c("vsearchMinQual", "vsearchMinovlen")])<=0)!=0){
    stop("Invalid 'vsearch...' arguments. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check positive integer fitnessHighConfidenceCount and fitnessDoubleHighConfidenceCount arguments
  if(sum(unlist(dimsum_meta[c("fitnessHighConfidenceCount", "fitnessDoubleHighConfidenceCount")])<0)!=0){
    stop("Invalid 'fitnessHighConfidenceCount' or 'fitnessDoubleHighConfidenceCount' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  #Check positive integer minimum read count arguments
  all_counts <- unlist(dimsum_meta[c("fitnessMinInputCountAll", "fitnessMinOutputCountAll", "fitnessMinInputCountAny", "fitnessMinOutputCountAny")])
  all_counts <- unlist(strsplit(all_counts, ",|:"))
  if(sum(!unique(unlist(strsplit(all_counts, ""))) %in% as.character(0:9))!=0){
    stop("Invalid 'fitness...Count...' arguments. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }else{
    #Parse minimum read count arguments
    for(i in c("fitnessMinInputCountAll", "fitnessMinOutputCountAll", "fitnessMinInputCountAny", "fitnessMinOutputCountAny")){
      dimsum_meta[[i]] <- suppressWarnings(dimsum__parse_minimum_read_count_arguments(dimsum_meta[[i]]))
    }
  }

  #Check strictly positive double splitChunkSize argument
  if(dimsum_meta[["splitChunkSize"]]<=0){
    stop("Invalid 'splitChunkSize' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check maxSubstitutions argument greater than 1
  if(dimsum_meta[["maxSubstitutions"]]<2){
    stop("Invalid 'maxSubstitutions' argument. Only integers greater than 1 allowed.", call. = FALSE)
  }

  #Check indels argument valid and set indelLengths
  dimsum_meta[["indelLengths"]] <- NA
  if(dimsum_meta[["indels"]]=="none"){
    #All indels discarded
    dimsum_meta[["indels"]] <- FALSE
  }else if(dimsum_meta[["indels"]]=="all"){
    #All indels retained (regardless of length)
    dimsum_meta[["indels"]] <- TRUE
  }else{
    #List of indel lengths specified
    dimsum_meta[["indelLengths"]] <- strtoi(unlist(strsplit(dimsum_meta[["indels"]], ",")))
    dimsum_meta[["indelLengths"]] <- dimsum_meta[["indelLengths"]][!is.na(dimsum_meta[["indelLengths"]])]
    if(length(dimsum_meta[["indelLengths"]])==0){
      stop("Invalid 'indels' argument. Only integers greater than 1 allowed.", call. = FALSE)
    }
    dimsum_meta[["indels"]] <- TRUE
  }

  #Check barcodeErrorRate argument
  if(dimsum_meta[["barcodeErrorRate"]]<0 | dimsum_meta[["barcodeErrorRate"]]>=1){
    stop("Invalid 'barcodeErrorRate' argument. Only positive doubles less than 1 allowed (zero inclusive).", call. = FALSE)
  }

  #Check vsearchMaxee argument
  if(dimsum_meta[["vsearchMaxee"]]<=0){
    stop("Invalid 'vsearchMaxee' argument. Only positive doubles allowed (zero exclusive).", call. = FALSE)
  }

  #Disable bayesianDoubleFitness option (TEMPORARY FIX: still in development)
  dimsum_meta[["bayesianDoubleFitness"]] <- FALSE

  #Check bayesianDoubleFitnessLamD argument
  if(dimsum_meta[["bayesianDoubleFitnessLamD"]]<=0){
    stop("Invalid 'bayesianDoubleFitnessLamD' argument. Only positive doubles allowed (zero exclusive).", call. = FALSE)
  }

  #Check startStage argument
  if(dimsum_meta[["startStage"]]<0){
    stop("Invalid 'startStage' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  #Check stopStage argument
  if(dimsum_meta[["stopStage"]]<=0){
    stop("Invalid 'stopStage' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
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
  if(!dimsum_meta[["sequenceType"]] %in% c("noncoding", "coding", "auto")){
    stop("Invalid 'sequenceType' argument. Only noncoding/coding/auto allowed.", call. = FALSE)
  }else if(dimsum_meta[["sequenceType"]]=="auto"){
    dimsum_meta[["sequenceType"]] <- dimsum__detect_sequence_type(dimsum_meta[["wildtypeSequence"]])
  }

  #Check mutagenesisType one of noncoding/coding/auto
  if(!dimsum_meta[["mutagenesisType"]] %in% c("random", "codon")){
    stop("Invalid 'mutagenesisType' argument. Only random/codon allowed.", call. = FALSE)
  }

  #Check retainedReplicates comma-separated list of integers or "all"
  if(dimsum_meta[["retainedReplicates"]] != "all"){
    if(sum(!strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]] %in% strsplit("0123456789", "")[[1]])!=0){
      stop("Invalid 'retainedReplicates' argument. Only comma-separated list of integers or 'all' allowed.", call. = FALSE)
    }
  }

  #Check if barcodeIdentityPath specified
  if(!is.null(dimsum_meta[["barcodeIdentityPath"]])){
    #Load barcode identity table
    if(!file.exists(dimsum_meta[["barcodeIdentityPath"]])){
      stop(paste0("Invalid '", "barcodeIdentityPath", "' argument (file not found)"), call. = FALSE)
    }
  }

  #Set fitnessNormalise to FALSE if fitnessErrorModel is FALSE
  if(!dimsum_meta[["fitnessErrorModel"]]){
    dimsum_meta[["fitnessNormalise"]] <- FALSE
  }

  #Set run start/end time
  dimsum_meta[["start_time"]] <- Sys.time()
  dimsum_meta[["end_time"]] <- "In progress"

  #Return
  return(dimsum_meta)
}
