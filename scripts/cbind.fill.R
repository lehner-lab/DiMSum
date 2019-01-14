
#' cbind.fill
#'
#' cbind a list of data.frames of same row names but unequal number.
#'
#' @param df_list a list of data.frames (required)
#'
#' @return a single data.frame where empty rows are filled with NAs
#' @export
cbind.fill <- function(df_list){
    nm <- lapply(df_list, as.matrix)
    n <- max(sapply(nm, nrow)) 
    temp_rownames <- unique(as.character(unlist(sapply(nm, rownames))))
    #Select longest rowname string if overlapping numerical ranges exist
    temp_rownames <- temp_rownames[order(nchar(temp_rownames), decreasing = T)]
    first_pos <- sapply(strsplit(temp_rownames, "-"), "[", 1)
    temp_rownames <- temp_rownames[!duplicated(first_pos)]
    first_pos <- as.numeric(first_pos[!duplicated(first_pos)])
    temp_rownames <- temp_rownames[order(first_pos, decreasing = F)]
    temp_df <- as.data.frame(do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))))
    rownames(temp_df) <- temp_rownames
    temp_df
  }
