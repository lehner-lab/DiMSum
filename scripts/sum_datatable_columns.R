
#' sum_datatable_columns
#'
#' Replace a subset of columns with a single column containing the row sums.
#'
#' @param dt Input data.table (required)
#' @param column_patterns character vector of column patterns to match (required)
#' @param suffix a character suffix for the appended column (default: "")
#'
#' @return a data.table where a subset of columns is replaced with a single column containing the row sums
#' @export
sum_datatable_columns <- function(
  dt, 
  column_patterns, 
  suffix=""){
  for(this_column_pattern in column_patterns){
    #Sum columns with given pattern
    temp_data <- apply(dt[,grep(this_column_pattern, colnames(dt)), with=FALSE], 1, sum)
    #Remove summed columns
    dt <- dt[,-grep(this_column_pattern, colnames(dt)), with=FALSE]
    #Append result column 
    dt[,paste0(this_column_pattern, suffix)] <- temp_data
  }
  #Return data.table
  return(dt)
}
