
#' dimsum_stage_vsearch
#'
#' Run VSEARCH on all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param vsearch_outpath VSEARCH output path (required)
#' @param report whether or not to generate VSEARCH summary plots (default: TRUE)
#' @param report_outpath VSEARCH report output path
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_vsearch <- function(
  dimsum_meta,
  vsearch_outpath,
  report = TRUE,
  report_outpath = NULL,
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

  #Create/overwrite vsearch directory (if executed)
  vsearch_outpath <- gsub("/$", "", vsearch_outpath)
  dimsum__create_dir(vsearch_outpath, execute = execute, message = "DiMSum STAGE 3 (WRAP): ALIGN PAIRED-END READS") 

  #Input files
  all_fastq <- file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], unique(c(dimsum_meta[['exp_design']][,"pair1"], dimsum_meta[['exp_design']][,"pair2"])))
  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = all_fastq,
    stage_number = this_stage,
    execute = execute)

  #Sample names
  sample_names <- paste0(
    dimsum_meta[["exp_design"]][,"sample_name"], '_e', 
    dimsum_meta[["exp_design"]][,"experiment"], '_s', 
    dimsum_meta[["exp_design"]][,"selection_id"], '_b', 
    dimsum_meta[["exp_design"]][,"biological_replicate"], '_t', 
    dimsum_meta[["exp_design"]][,"technical_replicate"], '_split', 
    dimsum_meta[["exp_design"]][,"split"], sep = "")

  #Additional vsearch options related to alignment length
  temp_options <- paste0(' -fastq_minovlen ', dimsum_meta[["vsearchMinovlen"]])
  if( dimsum_meta[["vsearchMinovlen"]] < 10 ){temp_options <- paste0(temp_options, " -fastq_maxdiffs ", dimsum_meta[["vsearchMinovlen"]])}
  # if( dimsum_meta[["vsearchMinovlen"]] < 16 ){temp_options <- paste0(temp_options, " -minhsp ", dimsum_meta[["vsearchMinovlen"]])}
  # if( dimsum_meta[["vsearchMinovlen"]] < 16 ){temp_options <- paste0(temp_options, " -band ", dimsum_meta[["vsearchMinovlen"]])}
  # if( dimsum_meta[["vsearchMinovlen"]] < 5 ){temp_options <- paste0(temp_options, " -hspw ", dimsum_meta[["vsearchMinovlen"]])}

  #Run VSEARCH on all fastq file pairs
  dimsum__status_message("Aligning paired-end FASTQ files with VSEARCH:\n")
  dimsum__status_message(paste0(all_fastq, "\n"))
  dimsum__status_message("Processing...\n")
  #Trans library mode?
  if(dimsum_meta[["transLibrary"]]){
    for(i in 1:length(sample_names)){dimsum__status_message(paste0("\t", paste0(unlist(dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')]), collapse = "\t"), "\n"))}
    if(execute){
      # Setup cluster
      clust <- parallel::makeCluster(dimsum_meta[['numCores']])
      # make variables available to each core's workspace
      parallel::clusterExport(clust, list("dimsum_meta","vsearch_outpath","sample_names","dimsum__concatenate_reads"), envir = environment())
      parallel::parSapply(clust,X = 1:length(sample_names), dimsum__vsearch_trans_library_helper)
      parallel::stopCluster(clust)
    }
  #Single-end mode?
  }else if(!dimsum_meta[["paired"]]){
    for(i in 1:length(sample_names)){dimsum__status_message(paste0("\t", paste0(unique(unlist(dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')])), collapse = "\t"), "\n"))}
    if(execute){
      # Setup cluster
      clust <- parallel::makeCluster(dimsum_meta[['numCores']])
      # make variables available to each core's workspace
      parallel::clusterExport(clust, list("dimsum_meta","vsearch_outpath","sample_names","dimsum__filter_single_end_reads"), envir = environment())
      parallel::parSapply(clust,X = 1:length(sample_names), dimsum__vsearch_single_end_library_helper)
      parallel::stopCluster(clust)
    }
  #Classic paired-end mode
  }else{
    for(i in 1:length(sample_names)){
      dimsum__status_message(paste0("\t", paste0(unlist(dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')]), collapse = "\t"), "\n"))
      #Check if this system command should be executed
      if(execute){
        temp_out <- system(paste0(
          "vsearch -fastq_mergepairs ",
          file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair1"]),
          " -reverse ",
          file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[["exp_design"]][i,"pair2"]),
          " -fastqout ",
          file.path(vsearch_outpath, paste0(sample_names[i], '.vsearch.prefilter')),
          " -quiet ",
          # " -fastq_minqual ",
          # as.character(dimsum_meta[["vsearchMinQual"]]),
          " -fastq_maxee ",
          as.character(dimsum_meta[["vsearchMaxee"]]),
          " -fastq_minlen ",
          as.character(dimsum_meta[["cutadaptMinLength"]]),
          temp_options,
          " -threads ",
          dimsum_meta[['numCores']],
          " --fastq_allowmergestagger ",
          " > ",
          file.path(vsearch_outpath, paste0(sample_names[i], '.vsearch.stdout')),
          " 2> ",
          file.path(vsearch_outpath, paste0(sample_names[i], '.report.prefilter'))))
      }
    }
    dimsum__status_message("Filtering aligned reads...\n")
    if(execute){
      #Input files
      input_files <- c(
        file.path(vsearch_outpath, paste0(sample_names, '.vsearch.prefilter')),
        file.path(vsearch_outpath, paste0(sample_names, '.report.prefilter')))
      #Check if all input files exist
      dimsum__check_files_exist(
        required_files = input_files,
        stage_number = this_stage,
        execute = execute)
      # Setup cluster
      clust <- parallel::makeCluster(dimsum_meta[['numCores']])
      # make variables available to each core's workspace
      parallel::clusterExport(clust, list("dimsum_meta","vsearch_outpath","sample_names","dimsum__filter_reads"), envir = environment())
      parallel::parSapply(clust,X = 1:length(sample_names), dimsum__filter_reads_helper)
      parallel::stopCluster(clust)
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Merged fastq filenames
  dimsum_meta_new[["exp_design"]][,"aligned_pair"] <- paste0(sample_names, ".vsearch")
  dimsum_meta_new[['exp_design']][,"aligned_pair_directory"] <- vsearch_outpath
  #Delete files when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      if(!is.null(dimsum_meta_new[["deleteIntermediateFiles"]])){suppressWarnings(temp_out <- file.remove(dimsum_meta_new[["deleteIntermediateFiles"]]))}
      if(!is.null(dimsum_meta_new[["touchIntermediateFiles"]])){suppressWarnings(temp_out <- file.create(dimsum_meta_new[["touchIntermediateFiles"]]))}
    }else{
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(
        dimsum_meta_new[["deleteIntermediateFiles"]], 
        file.path(vsearch_outpath, dir(vsearch_outpath, "*.vsearch$")),
        file.path(vsearch_outpath, dir(vsearch_outpath, "*.vsearch.prefilter$")))
    }
  }
  #Generate vsearch report
  if(report){
    tryCatch({
      dimsum_meta_new_report <- dimsum__vsearch_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
      }, error=function(e){
        dimsum__status_message("There were problems while running 'dimsum__vsearch_report'\n")
        dimsum_meta_new_report <- dimsum_meta_new
        })
    return(dimsum_meta_new_report)
  }
  return(dimsum_meta_new)
}
