
#dimsum_stage_unzip
#
# Unzip all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# fastq_outpath: FASTQ output path (required)
# execute: whether or not to execute the system command (default: TRUE)
# save_workspace: whether or not to save the current experiment metadata object (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_unzip <- function(
  dimsum_meta,
  fastq_outpath,
  execute = TRUE,
  save_workspace = TRUE
  ){
  #Create/overwrite unzip directory (if executed)
  fastq_outpath <- gsub("/$", "", fastq_outpath)
  create_dimsum_dir(fastq_outpath, execute = execute, message = "DiMSum STAGE 3: UNZIP")  
  #All fastq files gzipped?
  if(dimsum_meta[["gzipped"]]){
    message("Unzipping FASTQ files:")
    all_fastq <- file.path(dimsum_meta[["exp_design"]]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
    print(all_fastq)
    message("Processing...")
    for(f in all_fastq){
      message(paste0("\t", f))
      #Check if this system command should be executed
      if(execute){
        temp_out = system(paste0(
          "gunzip -c ", 
          f, 
          " > ", 
          file.path(fastq_outpath, gsub(".gz$", "", basename(f))),
          " 2> ",
          file.path(fastq_outpath, paste0(gsub(".gz$", "", basename(f)), '.stderr'))))
      }
    }
    #New experiment metadata
    dimsum_meta_new <- dimsum_meta
    #Update fastq metadata
    dimsum_meta_new[['exp_design']]$pair1 <- gsub(".gz$", "", dimsum_meta_new[["exp_design"]][,"pair1"])
    dimsum_meta_new[['exp_design']]$pair2 <- gsub(".gz$", "", dimsum_meta_new[["exp_design"]][,"pair2"])
    dimsum_meta_new[['exp_design']]$pair_directory <- fastq_outpath
    #Save workspace
    if(save_workspace){save_metadata(dimsum_meta_new)}
    return(dimsum_meta_new)
  }else{
    #Copy fastq files
    message("Skipping this stage (FASTQ files already unzipped)")
    #New experiment metadata
    dimsum_meta_new <- dimsum_meta
    #Save workspace
    if(save_workspace){save_metadata(dimsum_meta_new)}
    return(dimsum_meta_new)
  }
}

