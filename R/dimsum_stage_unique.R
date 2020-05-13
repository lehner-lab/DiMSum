
#' dimsum_stage_unique
#'
#' Run fastx_collapser on all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param unique_outpath fastx_collapser output path (required)
#' @param save_workspace whether or not to save the current experiment metadata object (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_unique <- function(
  dimsum_meta,
  unique_outpath,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 3
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #WRAP not run or this stage after stopStage
  if(!is.null(dimsum_meta[["countPath"]]) | dimsum_meta[["stopStage"]] < this_stage){
    return(dimsum_meta)
  }

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Create unique directory (if doesn't already exist)
  unique_outpath <- gsub("/$", "", unique_outpath)
  dimsum__create_dir(unique_outpath, execute = execute, message = "COUNT UNIQUE VARIANTS")  

  #Input files
  all_fasta <- file.path(dimsum_meta[["exp_design"]][,"aligned_pair_directory"], dimsum_meta[['exp_design']][,"aligned_pair"])
  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = all_fasta,
    stage_number = this_stage,
    execute = execute)

  #Add sample code (sample name plus experiment and replicate structure, but without split)
  dimsum_meta[["exp_design"]][,"sample_code"] <- sapply(strsplit(dimsum_meta[["exp_design"]][,"aligned_pair"], '.split'), '[', 1)

  #Run starcode on all aligned read pair fastq files
  dimsum__status_message("Counting unique aligned reads with starcode:\n")
  dimsum__status_message(paste0(all_fasta, "\n"))
  dimsum__status_message("Processing...\n")
  for(i in 1:dim(dimsum_meta[["exp_design"]])[1]){dimsum__status_message(paste0("\t", dimsum_meta[["exp_design"]][i,"aligned_pair"], "\n"))}
  #Check if this system command should be executed
  if(execute){
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","unique_outpath"), envir = environment())
    parallel::parSapply(clust,X = 1:length(unique(dimsum_meta[["exp_design"]][,"sample_code"])), dimsum__unique_helper)
    parallel::stopCluster(clust)
    #Run starcode to count unique variants
    for(this_sample_code in unique(dimsum_meta[["exp_design"]][,"sample_code"])){
      read_pairs <- dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"sample_code"]==this_sample_code,"aligned_pair"]
      output_file1 <- gsub("_split1.usearch$", ".usearch", read_pairs[1])
      output_file2 <- gsub("_split1.usearch$", ".usearch.unique", read_pairs[1])
      temp_out <- system(paste0(
        "starcode -d 0 -s -t ",
        dimsum_meta[['numCores']],
        " -i ",
        file.path(unique_outpath, output_file1),
        " -o ",
        file.path(unique_outpath, output_file2),
        " > ",
        file.path(unique_outpath, paste0(output_file2, '.stdout')),
        " 2> ",
        file.path(unique_outpath, paste0(output_file2, '.stderr'))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  dimsum_meta_new[["exp_design"]][,"aligned_pair_unique"] <- gsub("_split.*.usearch$", ".usearch.unique", dimsum_meta_new[["exp_design"]][,"aligned_pair"])
  dimsum_meta_new[['exp_design']][,"aligned_pair_unique_directory"] <- unique_outpath
  #Delete files when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      if(!is.null(dimsum_meta_new[["deleteIntermediateFiles"]])){suppressWarnings(temp_out <- file.remove(dimsum_meta_new[["deleteIntermediateFiles"]]))}
      if(!is.null(dimsum_meta_new[["touchIntermediateFiles"]])){suppressWarnings(temp_out <- file.create(dimsum_meta_new[["touchIntermediateFiles"]]))}
    }else{
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], 
        file.path(unique_outpath, dir(unique_outpath, "*.usearch$")),
        file.path(unique_outpath, dir(unique_outpath, "*.unique$")))
    }
  }
  return(dimsum_meta_new)
}

