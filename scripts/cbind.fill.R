
#cbind.fill
#
# cbind a list of data.frames of same row names but unequal number.
#
# df_list: a list of data.frames (required)
#
# Returns: a single data.frame where empty rows are filled with NAs.
#
cbind.fill <- function(df_list){
    nm <- lapply(df_list, as.matrix)
    n <- max(sapply(nm, nrow)) 
    temp_rownames <- unique(as.character(unlist(sapply(nm, rownames))))
    temp_df <- as.data.frame(do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))))
    rownames(temp_df) <- temp_rownames
    temp_df
  }
