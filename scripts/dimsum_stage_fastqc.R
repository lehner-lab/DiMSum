
#dimsum_stage_fastqc
#
# Run FASTQC on all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# fastqc_outpath: FASTQC output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_fastqc <- function(
  dimsum_meta,
  fastqc_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL
  ){
  #Create/overwrite FASTQC directory (if executed)
  fastqc_outpath <- gsub("/$", "", fastqc_outpath)
  create_dimsum_dir(fastqc_outpath, execute = execute, message = "DiMSum STAGE 2: FASTQC")  
  #Run FASTQC on all fastq files
  message("Running FASTQC on all files:")
  all_fastq <- file.path(dimsum_meta[['exp_design']]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
  print(all_fastq)
  message("Processing...")
  for(f in all_fastq){
    message(paste0("\t", f))
    #Check if this system command should be executed
    if(execute){
      temp_out = system(paste0(
        "fastqc -o ", 
        fastqc_outpath,
        " --extract ",
        " -t ",
        num_cores,
        " ",
        f,
        " > ", 
        file.path(fastqc_outpath, paste0(basename(f), '.stdout')),
        " 2> ",
        file.path(fastqc_outpath, paste0(basename(f), '.stderr'))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']]$pair1_fastqc <- gsub(dimsum_meta_new[["fastq_file_extension"]], '_fastqc/fastqc_data.txt', gsub('.gz', '', dimsum_meta_new[['exp_design']][,"pair1"]))
  dimsum_meta_new[['exp_design']]$pair2_fastqc <- gsub(dimsum_meta_new[["fastq_file_extension"]], '_fastqc/fastqc_data.txt', gsub('.gz', '', dimsum_meta_new[['exp_design']][,"pair2"]))
  dimsum_meta_new[['exp_design']]$fastqc_directory <- fastqc_outpath
  #Generate FASTQC report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_fastqc_report(dimsum_meta_new, report_outpath)
    return(dimsum_meta_new_report)
  }else{
    return(dimsum_meta_new)
  }
}

