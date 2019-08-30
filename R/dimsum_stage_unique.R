
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
  this_stage <- 5
  execute <- (dimsum_meta[["startStage"]] <= this_stage & (dimsum_meta[["stopStage"]] == 0 | dimsum_meta[["stopStage"]] >= this_stage))
  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}
  #Create unique directory (if doesn't already exist)
  unique_outpath <- gsub("/$", "", unique_outpath)
  dimsum__create_dir(unique_outpath, execute = execute, message = "DiMSum STAGE 5: COUNT UNIQUE VARIANTS")  
  #Run fastx_collapser on all aligned read pair fastq files
  message("Counting unique aligned reads with starcode:")
  all_fasta <- file.path(dimsum_meta[["exp_design"]][,"aligned_pair_directory"], dimsum_meta[['exp_design']][,"aligned_pair"])
  print(all_fasta)
  message("Processing...")
  for(i in 1:dim(dimsum_meta[["exp_design"]])[1]){message(paste0("\t", dimsum_meta[["exp_design"]][i,"aligned_pair"]))}
  #Check if this system command should be executed
  if(execute){
    dimsum_stage_unique_helper <- function(
      i
      ){
      this_sample_name <- unique(dimsum_meta[["exp_design"]][,"sample_name"])[i]
      read_pairs <- dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"sample_name"]==this_sample_name,"aligned_pair"]
      #Concatenate FASTQ files
      output_file1 <- gsub("_split1.usearch$", ".usearch", read_pairs[1])
      temp_out <- system(paste0(
        "cat ",
        paste(file.path(dimsum_meta[["exp_design"]][,"aligned_pair_directory"][1], read_pairs), collapse = " "),
        " > ",
        file.path(unique_outpath, output_file1)))
    }
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","unique_outpath"), envir = environment())
    parallel::parSapply(clust,X = 1:length(unique(dimsum_meta[["exp_design"]][,"sample_name"])), dimsum_stage_unique_helper)
    parallel::stopCluster(clust)
    #Run starcode to count unique variants
    for(this_sample_name in unique(dimsum_meta[["exp_design"]][,"sample_name"])){
      read_pairs <- dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"sample_name"]==this_sample_name,"aligned_pair"]
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
      temp_out <- mapply(system, dimsum_meta_new[["deleteIntermediateFiles"]], MoreArgs = list(ignore.stdout = T, ignore.stderr = T))
    }else{
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], 
        paste0("rm ", file.path(unique_outpath, "*.usearch")),
        paste0("rm ", file.path(unique_outpath, "*.unique")))
    }
  }
  return(dimsum_meta_new)
}

