
#' dimsum__parse_minimum_read_count_arguments
#'
#' Parse minimum input count arguments.
#'
#' @param input_arg a string (required)
#'
#' @return a single integer or list of named integers
#' @export
dimsum__parse_minimum_read_count_arguments <- function(
  input_arg
  ){
  #A single integer supplied
  if(!is.na(as.integer(input_arg))){
    return(as.integer(input_arg))
  }

  #Multiple integers supplied
  temp_list <- sapply(as.list(unlist(strsplit(input_arg, ","))), strsplit, ":")
  names(temp_list) <- sapply(temp_list, '[', 1)
  #List of list of 2 integers
  if(sum(is.na(as.integer(unlist(temp_list))))!=0 | sum(sapply(temp_list, length)!=2)!=0 | sum(duplicated(names(temp_list)))!=0){
    stop("Invalid 'fitness...Count...' arguments. Multiple thresholds should be specified in the form 'a:b,c:d' where a,b,c and d are positive integers.", call. = FALSE)
  }
  #Set names
  names(temp_list) <- sapply(temp_list, '[', 1)
  #Return list
  return(lapply(lapply(temp_list, '[', 2), as.integer))
}

