
#' dimsum__cbind_fill
#'
#' cbind a list of data.frames of same row names but unequal number.
#'
#' @param df_list a list of data.frames (required)
#'
#' @return a single data.frame where empty rows are filled with NAs
#' @export
dimsum__cbind_fill <- function(df_list){
    nm <- lapply(lapply(df_list, dimsum__derange_df), as.matrix)
    n <- max(sapply(nm, nrow)) 
    temp_df <- as.data.frame(do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))))
    rownames(temp_df) <- 1:dim(temp_df)[1]
    temp_df <- dimsum__range_df(temp_df)
    colnames(temp_df) <- names(df_list)
    temp_df
  }
