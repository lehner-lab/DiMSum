
#' dimsum_stage_unique
#'
#' Run fastx_collapser on all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param unique_outpath fastx_collapser output path (required)
#' @param execute whether or not to execute the system command (default: TRUE)
#' @param save_workspace whether or not to save the current experiment metadata object (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_unique <- function(
  dimsum_meta,
  unique_outpath,
  execute = TRUE,
  save_workspace = TRUE
  ){
  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}
  #Create unique directory (if doesn't already exist)
  unique_outpath <- gsub("/$", "", unique_outpath)
  dimsum__create_dir(unique_outpath, execute = execute, message = "DiMSum STAGE 5: COUNT UNIQUE VARIANTS")  
  #Run fastx_collapser on all aligned read pair fastq files
  message("Getting unique aligned read counts with fastx_collapser:")
  all_fasta <- file.path(dimsum_meta[["exp_design"]][,"aligned_pair_directory"], dimsum_meta[['exp_design']][,"aligned_pair"])
  print(all_fasta)
  message("Processing...")
  for(i in 1:dim(dimsum_meta[["exp_design"]])[1]){message(paste0("\t", dimsum_meta[["exp_design"]][i,"aligned_pair"]))}
  #Check if this system command should be executed
  if(execute){
    dimsum_stage_unique_helper <- function(
      i
      ){
      read_pair <- dimsum_meta[["exp_design"]][i,"aligned_pair"]
      temp_out = system(paste0(
        "fastx_collapser -Q33 -i ",
        file.path(dimsum_meta[["exp_design"]][,"aligned_pair_directory"][1], read_pair),
        " -o ",
        file.path(unique_outpath, paste0(read_pair, '.unique')),
        " > ",
        file.path(unique_outpath, paste0(read_pair, '.unique.stdout')),
        " 2> ",
        file.path(unique_outpath, paste0(read_pair, '.unique.stderr'))))
    }
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","unique_outpath"), envir = environment())
    parallel::parSapply(clust,X = 1:nrow(dimsum_meta[["exp_design"]]), dimsum_stage_unique_helper)
    parallel::stopCluster(clust)
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Unique fasta filenames
  dimsum_meta_new[["exp_design"]][,"aligned_pair_unique"] <- paste0(dimsum_meta_new[["exp_design"]][,"aligned_pair"], ".unique")
  dimsum_meta_new[['exp_design']][,"aligned_pair_unique_directory"] <- unique_outpath
  return(dimsum_meta_new)
}

