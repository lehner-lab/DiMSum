
#' dimsum__check_variants
#'
#' Check whether minimum variant requirements satisfied.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__check_variants <- function(
  dimsum_meta,
  input_dt
  ){

  #Determine replicates to retain
  all_reps <- unique(dimsum_meta[["exp_design"]][,"experiment"])
  if(dimsum_meta[["retainedReplicates"]]!="all"){
    all_reps <- unique(as.integer(strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]]))
  }

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="|")

  #Check if WT sequence and at least 2 substitutions variants exist
  if(input_dt[WT==T,.N]==0){
    
    if( length(grep(paste0("count_e(", all_reps_str, ")_s[01]"), names(input_dt))) != 0 ){

	    input_dt[, all_reads := rowSums(.SD > 0) == (2*length(all_reps)),,.SDcols = grep(paste0("count_e(", all_reps_str, ")_s[01]"), names(input_dt))]
	    input_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("count_e(", all_reps_str, ")_s0"), names(input_dt))]
	    dimsum__status_message(paste0("WT variant has zero count in at least one input/output replicate. Did you mean to specify one of the following?\n"))
	    if(dimsum_meta[["sequenceType"]]=="coding" & dimsum_meta[["mixedSubstitutions"]]){
	      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(aa_seq, all_reads, mean_count)])
	    }else if(dimsum_meta[["sequenceType"]]=="coding"){
	      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), aa_seq, all_reads, mean_count)])
	    }else{
	      print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
	    }
	    stop(paste0("Cannot proceed with fitness estimation: WT variant not found"), call. = FALSE)

    }else{

	    input_dt[, all_reads := rowSums(.SD > 0) == length(grep(paste0("e(", all_reps_str, ")_s[01].*_count$"), names(input_dt))),,.SDcols = grep(paste0("e(", all_reps_str, ")_s[01].*_count$"), names(input_dt))]
	    input_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("e(", all_reps_str, ")_s[0].*_count$"), names(input_dt))]
	    dimsum__status_message(paste0("WT variant not found. Did you mean to specify one of the following?\n"))
	    print(input_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
	    stop(paste0("Cannot proceed with fitness estimation: WT variant not found"), call. = FALSE)  

    }

  }else if(input_dt[is.na(WT),.N]<2){
    stop(paste0("Cannot proceed with fitness estimation: insufficient substitution variants"), call. = FALSE)
  }

}
