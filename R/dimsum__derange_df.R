
#' dimsum__derange_df
#'
#' Expand a data.frame with ranged row.names.
#'
#' @param input_df input data.frame (required)
#'
#' @return a single deranged data.frame 
#' @export
dimsum__derange_df <- function(input_df){
    all_ind <- rownames(input_df)
    max_ind <- max(as.numeric(unlist(strsplit(all_ind[length(all_ind)], "-"))))
    output_list <- as.list(rep(NA, max_ind))
    max_ind_all <- sapply(lapply(lapply(strsplit(all_ind, "-"), unlist), as.numeric), min)
    output_list[max_ind_all] <- input_df[,1]
    for(i in 1:length(output_list)){
        if(is.na(output_list[[i]])){
            if(i==1){
                output_list[[i]] <- NA
            }else{
                output_list[[i]] <- output_list[[i-1]]
            }
        }
    }
    output_df <- do.call("rbind", output_list)
    rownames(output_df) <- 1:length(output_list)
    colnames(output_df) <- colnames(input_df)
    return(output_df)
}
