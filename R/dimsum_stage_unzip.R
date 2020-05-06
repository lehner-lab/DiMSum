
#' dimsum_stage_unzip
#'
#' Unzip all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param fastq_outpath FASTQ output path (required)
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_unzip <- function(
  dimsum_meta,
  fastq_outpath,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 1
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #WRAP not run
  if(!is.null(dimsum_meta[["countPath"]])){
    return(dimsum_meta)
  }

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Create/overwrite unzip directory (if executed)
  fastq_outpath <- gsub("/$", "", fastq_outpath)
  dimsum__create_dir(fastq_outpath, execute = execute, message = "UNZIP FASTQ FILES")  

  #All fastq files gzipped?
  if(dimsum_meta[["gzipped"]]){
    #Input files
    all_fastq <- file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], unique(c(dimsum_meta[['exp_design']][,"pair1"], dimsum_meta[['exp_design']][,"pair2"])))
    #Check if all input files exist
    dimsum__check_files_exist(
      required_files = all_fastq,
      stage_number = this_stage,
      execute = execute)

    dimsum__status_message("Unzipping FASTQ files:\n")
    dimsum__status_message(paste0(all_fastq, "\n"))
    dimsum__status_message("Processing...\n")
    dimsum__status_message(paste0("\t", basename(all_fastq), "\n"))
    #Check if this system command should be executed
    if(execute){
      # Setup cluster
      clust <- parallel::makeCluster(dimsum_meta[['numCores']])
      # make variables available to each core's workspace
      parallel::clusterExport(clust, list("all_fastq","fastq_outpath"), envir = environment())
      parallel::parSapply(clust,X = 1:length(all_fastq), dimsum__unzip_helper)
      parallel::stopCluster(clust)
    }
    #New experiment metadata
    dimsum_meta_new <- dimsum_meta
    #Update fastq metadata
    dimsum_meta_new[['exp_design']][,"pair1"] <- gsub(".gz$", "", dimsum_meta_new[["exp_design"]][,"pair1"])
    dimsum_meta_new[['exp_design']][,"pair2"] <- gsub(".gz$", "", dimsum_meta_new[["exp_design"]][,"pair2"])
    dimsum_meta_new[['exp_design']][,"pair_directory"] <- fastq_outpath
    #Delete files when last stage complete
    if(!dimsum_meta_new[["retainIntermediateFiles"]]){
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], 
        file.path(dimsum_meta_new[['exp_design']][,"pair_directory"], c(dimsum_meta_new[['exp_design']][,"pair1"], dimsum_meta_new[['exp_design']][,"pair2"])))
    }
    return(dimsum_meta_new)
  }
  #Copy fastq files
  dimsum__status_message("Skipping this stage (FASTQ files already unzipped)\n")
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  return(dimsum_meta_new)
}

