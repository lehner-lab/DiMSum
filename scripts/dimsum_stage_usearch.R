
#dimsum_stage_usearch
#
# Run USEARCH on all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# usearch_outpath: USEARCH output path (required)
# execute: whether or not to execute the system command (default: TRUE)
# save_workspace: whether or not to save the current experiment metadata object (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_usearch <- function(
  dimsum_meta,
  usearch_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){
  #Create/overwrite usearch directory (if executed)
  usearch_outpath <- gsub("/$", "", usearch_outpath)
  create_dimsum_dir(usearch_outpath, execute = execute, message = "DiMSum STAGE 6: USEARCH")  
  #Sample names
  sample_names = paste0(
    dimsum_meta[["exp_design"]][,"sample_name"], '_e', 
    dimsum_meta[["exp_design"]][,"experiment"], '_s', 
    dimsum_meta[["exp_design"]][,"selection_id"], '_b', 
    dimsum_meta[["exp_design"]][,"biological_replicate"], '_t', 
    dimsum_meta[["exp_design"]][,"technical_replicate"], '_split', 
    dimsum_meta[["exp_design"]][,"split"], sep = "")
  #Additional usearch options related to alignment length
  temp_options = paste0(' -fastq_minovlen ', dimsum_meta[["usearchMinovlen"]])
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options = paste0(temp_options, " -xdrop_nw ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options = paste0(temp_options, " -minhsp ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options = paste0(temp_options, " -band ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 5 ){temp_options = paste0(temp_options, " -hspw ", dimsum_meta[["usearchMinovlen"]])}
  #Run USEARCH on all fastq file pairs
  message("Aligning paired-end FASTQ files with USEARCH:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]][,"pair_directory"], c(dimsum_meta[['exp_design']][,"pair1"], dimsum_meta[['exp_design']][,"pair2"]))
  print(all_fastq)
  message("Processing...")
  for(i in 1:length(sample_names)){
    #TODO: usearch binary path specifiable on commandline?
    #TODO: only run if usearch arguments specified
    message(paste0("\t", dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')]))
    #Check if this system command should be executed
    if(execute){
      if(dimsum_meta[["transLibrary"]]){
        temp_out = fastq_manualalign(
          input_FASTQ1 = file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair1"]),
          input_FASTQ2 = file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair2"]),
          output_FASTQ = file.path(usearch_outpath, paste0(sample_names[i], '.usearch')),
          output_REPORT = file.path(usearch_outpath, paste0(sample_names[i], '.report')),
          num_nuc = 0,
          min_qual = dimsum_meta[["usearchMinQual"]],
          max_ee = dimsum_meta[["usearchMaxee"]],
          min_len = dimsum_meta[["usearchMinlen"]],
          concatentate_reads = TRUE)              
      }
      else if(dimsum_meta[["usearchAttemptExactMinovlen"]]){
        temp_out = system(paste0(
          "fastq_manualalign.py -i1 ",
          file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair1"]),
          " -i2 ",
          file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair2"]),
          " -o ",
          file.path(usearch_outpath, paste0(sample_names[i], '.usearch')),
          " -r ",
          file.path(usearch_outpath, paste0(sample_names[i], '.report')),
          " --minqual ",
          as.character(dimsum_meta[["usearchMinQual"]]),
          " --maxee ",
          as.character(dimsum_meta[["usearchMaxee"]]),
          " --numNuc ",
          as.character(dimsum_meta[["usearchMinovlen"]]),
          " > ",
          file.path(usearch_outpath, paste0(sample_names[i], '.usearch.stdout')),
          " 2> ",
          file.path(usearch_outpath, paste0(sample_names[i], '.usearch.stderr'))))        
      }
      else{
        temp_out = system(paste0(
          "usearch -fastq_mergepairs ",
          file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair1"]),
          " -reverse ",
          file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair2"]),
          " -fastqout ",
          file.path(usearch_outpath, paste0(sample_names[i], '.usearch')),
          " -report ",
          file.path(usearch_outpath, paste0(sample_names[i], '.report')),
          " -fastq_minqual ",
          as.character(dimsum_meta[["usearchMinQual"]]),
          " -fastq_merge_maxee ",
          as.character(dimsum_meta[["usearchMaxee"]]),
          " -fastq_minlen ",
          as.character(dimsum_meta[["usearchMinlen"]]),
          temp_options,
          " -threads ",
          dimsum_meta[['num_cores']],
          " > ",
          file.path(usearch_outpath, paste0(sample_names[i], '.usearch.stdout')),
          " 2> ",
          file.path(usearch_outpath, paste0(sample_names[i], '.usearch.stderr'))))
      }
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Merged fastq filenames
  dimsum_meta_new[["exp_design"]][,"aligned_pair"] <- paste0(sample_names, ".usearch")
  dimsum_meta_new[['exp_design']][,"aligned_pair_directory"] <- usearch_outpath
  #Generate usearch report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_usearch_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
    #Save workspace
    if(save_workspace){save_metadata(dimsum_meta_new_report)}
    return(dimsum_meta_new_report)
  }
  #Save workspace
  if(save_workspace){save_metadata(dimsum_meta_new)}
  return(dimsum_meta_new)
}
