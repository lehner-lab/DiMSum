
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
  if(save_workspace){save_metadata(dimsum_meta = dimsum_meta, n = 2)}
  #Create unique directory (if doesn't already exist)
  unique_outpath <- gsub("/$", "", unique_outpath)
  create_dimsum_dir(unique_outpath, execute = execute, message = "DiMSum STAGE 5: COUNT UNIQUE VARIANTS")  
  #Run fastx_collapser on all aligned read pair fastq files
  message("Getting unique aligned read counts with fastx_collapser:")
  all_fasta <- file.path(dimsum_meta[["exp_design"]][,"aligned_pair_directory"], dimsum_meta[['exp_design']][,"aligned_pair"])
  print(all_fasta)
  message("Processing...")
  for(read_pair in dimsum_meta[["exp_design"]][,"aligned_pair"]){
    #TODO: usearch binary path specifiable on commandline?
    #TODO: only run if usearch arguments specified
    message(paste0("\t", read_pair))
    #Check if this system command should be executed
    if(execute){
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
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Unique fasta filenames
  dimsum_meta_new[["exp_design"]][,"aligned_pair_unique"] <- paste0(dimsum_meta_new[["exp_design"]][,"aligned_pair"], ".unique")
  dimsum_meta_new[['exp_design']][,"aligned_pair_unique_directory"] <- unique_outpath
  return(dimsum_meta_new)
}

