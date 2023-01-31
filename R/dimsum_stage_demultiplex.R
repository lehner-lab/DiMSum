
#' dimsum_stage_demultiplex
#'
#' Run demultiplexing on all fastq files using cutadapt.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param demultiplex_outpath demultiplex output path (required)
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_demultiplex <- function(
  dimsum_meta,
  demultiplex_outpath,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 0
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #WRAP not run
  if(!is.null(dimsum_meta[["countPath"]])){
    return(dimsum_meta)
  }

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Create/overwrite demultiplex directory (if executed)
  demultiplex_outpath <- gsub("/$", "", demultiplex_outpath)
  dimsum__create_dir(demultiplex_outpath, execute = execute, message = "DiMSum STAGE 0 (WRAP): DEMULTIPLEX READS")  

  #Demultiplex parameters specified?
  if( !'barcode_design' %in% names(dimsum_meta) ){
    dimsum__status_message("Skipping this stage (assuming all fastq files already demultiplexed)\n")
    return(dimsum_meta)
  }

  #Input files
  fastq_pair_list <- unique(dimsum_meta[['barcode_design']][,c('pair1', 'pair2')])
  rownames(fastq_pair_list) = 1:dim(fastq_pair_list)[1]
  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], unlist(fastq_pair_list)),
    stage_number = this_stage,
    execute = execute)
  
  #Reformat barcode files for cutadapt (and copy and rename FASTQ files if necessary)
  #Check if this system command should be executed
  if(execute){
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","fastq_pair_list","demultiplex_outpath"), envir = environment())
    parallel::parSapply(clust,X = 1:nrow(fastq_pair_list), dimsum__demultiplex_cp_helper)
    parallel::stopCluster(clust)
  }
  #Update names in list if file extension incompatible with cutadapt (i.e. NOT ".fastq" or ".fastq.gz")
  if(dimsum_meta[["fastqFileExtension"]]!=".fastq"){
    for(pair_name in rownames(fastq_pair_list)){
      #New FASTQ file names
      new_fastq_name1 <- gsub(paste0(dimsum_meta[["fastqFileExtension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), ".fastq.gz", fastq_pair_list[pair_name,][1])
      new_fastq_name2 <- gsub(paste0(dimsum_meta[["fastqFileExtension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), ".fastq.gz", fastq_pair_list[pair_name,][2])
      #Update names in list
      fastq_pair_list[pair_name,][1] <- new_fastq_name1
      fastq_pair_list[pair_name,][2] <- new_fastq_name2
    }
    #Update FASTQ file directory in list
    dimsum_meta[["exp_design"]][,"pair_directory"] <- demultiplex_outpath
  }
  #Demultiplex FASTQ files
  dimsum__status_message("Demultiplexing FASTQ files with cutadapt:\n")
  all_fastq <- file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], unlist(fastq_pair_list))
  dimsum__status_message(paste0(unique(all_fastq), "\n"))
  dimsum__status_message("Processing...\n")
  for(i in 1:dim(fastq_pair_list)[1]){dimsum__status_message(paste0("\t", unique(unlist(fastq_pair_list[i,])), "\n"))}
  #Check if this system command should be executed
  if(execute){
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","demultiplex_outpath","fastq_pair_list"), envir = environment())
    parallel::parSapply(clust,X = 1:nrow(fastq_pair_list), dimsum__demultiplex_helper)
    parallel::stopCluster(clust)
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Delete files when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], file.path(demultiplex_outpath, dir(demultiplex_outpath, "*.fastq.gz$")))
    if(dimsum_meta[["fastqFileExtension"]]!=".fastq"){
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], unique(all_fastq))
    }
  }
  #Update fastq metadata
  dimsum_meta_new[['exp_design']][,"pair_directory"] <- demultiplex_outpath
  dimsum_meta_new[["fastqFileExtension"]] <- ".fastq"
  return(dimsum_meta_new)
}

