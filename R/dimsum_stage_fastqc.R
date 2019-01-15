
#' dimsum_stage_fastqc
#'
#' Run FASTQC on all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param fastqc_outpath FASTQC output path (required)
#' @param execute whether or not to execute the system command (default: TRUE)
#' @param report whether or not to generate FASTQC summary plots (default: TRUE)
#' @param report_outpath FASTQC report output path
#' @param save_workspace whether or not to save the current experiment metadata object (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_fastqc <- function(
  dimsum_meta,
  fastqc_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){
  #Create/overwrite FASTQC directory (if executed)
  fastqc_outpath <- gsub("/$", "", fastqc_outpath)
  create_dimsum_dir(fastqc_outpath, execute = execute, message = "DiMSum STAGE 2: ASSESS READ QUALITY")  
  #Run FASTQC on all fastq files
  message("Running FASTQC on all files:")
  all_fastq <- file.path(dimsum_meta[['exp_design']][,"pair_directory"], c(dimsum_meta[['exp_design']][,'pair1'], dimsum_meta[['exp_design']][,'pair2']))
  print(all_fastq)
  message("Processing...")
  message(paste0("\t", all_fastq, "\n"))
  #Check if this system command should be executed
  if(execute){
    temp_out = system(paste0(
      "fastqc -o ", 
      fastqc_outpath,
      " --extract ",
      " -t ",
      dimsum_meta[['num_cores']],
      " ",
      paste(all_fastq, collapse = " "),
      " > ", 
      file.path(fastqc_outpath, paste0('fastqc', '.stdout')),
      " 2> ",
      file.path(fastqc_outpath, paste0('fastqc', '.stderr'))))
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']][,"pair1_fastqc"] <- gsub(dimsum_meta_new[["fastq_file_extension"]], '_fastqc/fastqc_data.txt', gsub('.gz', '', dimsum_meta_new[['exp_design']][,"pair1"]))
  dimsum_meta_new[['exp_design']][,"pair2_fastqc"] <- gsub(dimsum_meta_new[["fastq_file_extension"]], '_fastqc/fastqc_data.txt', gsub('.gz', '', dimsum_meta_new[['exp_design']][,"pair2"]))
  dimsum_meta_new[['exp_design']][,"fastqc_directory"] <- fastqc_outpath
  #Generate FASTQC report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_fastqc_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
    #Save workspace
    if(save_workspace){save_metadata(dimsum_meta_new_report)}
    return(dimsum_meta_new_report)
  }
  #Save workspace
  if(save_workspace){save_metadata(dimsum_meta_new)}
  return(dimsum_meta_new)
}

