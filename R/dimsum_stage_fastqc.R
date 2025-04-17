
#' dimsum_stage_fastqc
#'
#' Run FASTQC on all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param fastqc_outpath FASTQC output path (required)
#' @param report whether or not to generate FASTQC summary plots (default: TRUE)
#' @param report_outpath FASTQC report output path
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_fastqc <- function(
  dimsum_meta,
  fastqc_outpath,
  report = TRUE,
  report_outpath = NULL,
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

  #Create/overwrite FASTQC directory (if executed)
  fastqc_outpath <- gsub("/$", "", fastqc_outpath)
  dimsum__create_dir(fastqc_outpath, execute = execute, message = "DiMSum STAGE 1 (WRAP): ASSESS READ QUALITY")  

  #Input files
  all_fastq <- unique(c(
    file.path(dimsum_meta[['exp_design']][,"pair_directory"], dimsum_meta[['exp_design']][,'pair1']),
    file.path(dimsum_meta[['exp_design']][,"pair_directory"], dimsum_meta[['exp_design']][,'pair2'])))

  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = all_fastq,
    stage_number = this_stage,
    execute = execute)

  #Run FASTQC on all fastq files
  dimsum__status_message("Running FASTQC on all files:\n")
  dimsum__status_message(paste0(all_fastq, "\n"))
  dimsum__status_message("Processing...\n")
  dimsum__status_message(paste0("\t", basename(all_fastq), "\n"))
  #Check if this system command should be executed
  if(execute){
    temp_out <- system(paste0(
      "fastqc -o ", 
      fastqc_outpath,
      " --extract ",
      " -t ",
      dimsum_meta[['numCores']],
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
  dimsum_meta_new[['exp_design']][,"pair1_fastqc"] <- gsub(dimsum_meta_new[["fastqFileExtension"]], '_fastqc/fastqc_data.txt', gsub('.gz', '', dimsum_meta_new[['exp_design']][,"pair1"]))
  dimsum_meta_new[['exp_design']][,"pair2_fastqc"] <- gsub(dimsum_meta_new[["fastqFileExtension"]], '_fastqc/fastqc_data.txt', gsub('.gz', '', dimsum_meta_new[['exp_design']][,"pair2"]))
  dimsum_meta_new[['exp_design']][,"fastqc_directory"] <- fastqc_outpath
  #Generate FASTQC report
  if(report){
    tryCatch({
      dimsum_meta_new_report <- dimsum__fastqc_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
      }, error=function(e){
        dimsum__status_message("There were problems while running 'dimsum__fastqc_report'\n")
        dimsum_meta_new_report <- dimsum_meta_new
        })
    return(dimsum_meta_new_report)
  }
  return(dimsum_meta_new)
}

