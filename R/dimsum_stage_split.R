
#' dimsum_stage_split
#'
#' Split all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param split_outpath split FASTQ output path (required)
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_split <- function(
  dimsum_meta,
  split_outpath,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 1
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #WRAP not run or this stage after stopStage
  if(!is.null(dimsum_meta[["countPath"]]) | dimsum_meta[["stopStage"]] < this_stage){
    return(dimsum_meta)
  }

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Create/overwrite split directory (if executed)
  split_outpath <- gsub("/$", "", split_outpath)
  dimsum__create_dir(split_outpath, execute = execute, message = "SPLIT FASTQ FILES")  

  #Input files
  fastq_pair_list <- unique(dimsum_meta[['exp_design']][,c('pair1', 'pair2')])
  rownames(fastq_pair_list) <- 1:dim(fastq_pair_list)[1]
  all_fastq <- file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], unique(c(dimsum_meta[['exp_design']][,"pair1"], dimsum_meta[['exp_design']][,"pair2"])))
  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = all_fastq,
    stage_number = this_stage,
    execute = execute)

  #Split FASTQ files
  dimsum__status_message("Splitting FASTQ files:\n")
  dimsum__status_message(paste0(all_fastq, "\n"))
  dimsum__status_message("Processing...\n")
  dimsum__status_message(paste0("\t", basename(all_fastq), "\n"))
  #Check if this system command should be executed
  if(execute){
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","fastq_pair_list","split_outpath","dimsum__fastq_splitter","dimsum__writeFastq"), envir = environment())
    parallel::parSapply(clust,X = 1:nrow(fastq_pair_list), dimsum__split_helper)
    parallel::stopCluster(clust)
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  all_fastq <- file.path(dimsum_meta_new[["exp_design"]][,"pair_directory"], c(dimsum_meta_new[['exp_design']][,"pair1"]))
  split_list <- list()
  for(f in all_fastq){
    num_splits <- length(list.files(split_outpath, pattern = basename(f)))
    split_list <- append(split_list, num_splits)
  }
  dimsum_meta_new[["exp_design"]] <- dimsum_meta_new[["exp_design"]][rep(1:length(all_fastq), times = unlist(split_list)),]
  temp_rownames <- rownames(dimsum_meta_new[["exp_design"]])
  temp_suffix <- rep('.split1', dim(dimsum_meta_new[["exp_design"]])[1])
  temp_suffix[grepl('\\.', temp_rownames)] <- paste0('.split', as.integer(sapply(strsplit(temp_rownames[grepl('\\.', temp_rownames)], '\\.'), '[', 2))+1)
  dimsum_meta_new[["exp_design"]][,"pair1"] <- paste0(dimsum_meta_new[["exp_design"]][,"pair1"], temp_suffix, '.fastq')
  dimsum_meta_new[["exp_design"]][,"pair2"] <- paste0(dimsum_meta_new[["exp_design"]][,"pair2"], temp_suffix, '.fastq')
  dimsum_meta_new[["exp_design"]][,"split"] <- as.integer(gsub(".split", "", temp_suffix))
  dimsum_meta_new[['exp_design']][,"pair_directory"] <- split_outpath
  #Delete file contents when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      if(!is.null(dimsum_meta_new[["deleteIntermediateFiles"]])){suppressWarnings(temp_out <- file.remove(dimsum_meta_new[["deleteIntermediateFiles"]]))}
      if(!is.null(dimsum_meta_new[["touchIntermediateFiles"]])){suppressWarnings(temp_out <- file.create(dimsum_meta_new[["touchIntermediateFiles"]]))}
    }else{
      dimsum_meta_new[["touchIntermediateFiles"]] <- file.path(dimsum_meta_new[['exp_design']][,"pair_directory"], c(dimsum_meta_new[['exp_design']][,"pair1"], dimsum_meta_new[['exp_design']][,"pair2"]))
    }
  }
  return(dimsum_meta_new)
}

