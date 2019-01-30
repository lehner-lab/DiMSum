
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

  #Check WT sequence is valid (ACGT characters only)
  all_characters <- unique(unlist(strsplit(dimsum_meta[["wildtypeSequence"]], "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T"))!=0){
    stop("Invalid wild-type nucleotide sequence. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  #Check strictly positive integer usearch... arguments
  if(sum(unlist(dimsum_meta[c("usearchMinQual", "usearchMinlen", "usearchMinovlen")])<=0)!=0){
    stop("Invalid 'usearch...' arguments. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }

  #Check usearchMaxee argument
  if(dimsum_meta[["usearchMaxee"]]<=0){
    stop("Invalid 'usearchMaxee' argument. Only positive doubles allowed (zero exclusive).", call. = FALSE)
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

  #Return
  return(dimsum_meta)
}
