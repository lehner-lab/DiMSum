
#' dimsum__range_df
#'
#' Collapse a data.frame with integer rownames to ranged row.names.
#'
#' @param input_df input data.frame (required)
#'
#' @return a single ranged data.frame 
#' @export
dimsum__range_df <- function(input_df){
    output_df <- input_df[!duplicated(input_df),]
    rn <- as.numeric(rownames(output_df))
    rn_list <- as.list(rn)
    for(i in 1:length(rn)){
        if(i!=1){
            if(rn_list[[i]]>(rn_list[[i-1]]+1)){
                rn_list[[i-1]] <- c(rn_list[[i-1]], rn_list[[i]]-1)
            }
        }
    }
    ind_range <- sapply(lapply(lapply(rn_list, range), unique), paste, collapse="-")
    #Correct for terminal duplicated rows
    if(max(as.numeric(unlist(strsplit(ind_range, "-"))))!=dim(input_df)[1]){
        final_range <- unique(range(c(as.numeric(unlist(strsplit(ind_range[length(ind_range)], "-"))), dim(input_df)[1])))
        ind_range[length(ind_range)] <- paste(final_range, collapse = "-")
    }
    rownames(output_df) <- ind_range
    return(output_df)
}
