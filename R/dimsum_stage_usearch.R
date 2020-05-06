
#' dimsum_stage_usearch
#'
#' Run USEARCH on all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param usearch_outpath USEARCH output path (required)
#' @param report whether or not to generate USEARCH summary plots (default: TRUE)
#' @param report_outpath USEARCH report output path
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_usearch <- function(
  dimsum_meta,
  usearch_outpath,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 3
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #WRAP not run
  if(!is.null(dimsum_meta[["countPath"]])){
    return(dimsum_meta)
  }

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Create/overwrite usearch directory (if executed)
  usearch_outpath <- gsub("/$", "", usearch_outpath)
  dimsum__create_dir(usearch_outpath, execute = execute, message = "DiMSum STAGE 3: ALIGN PAIRED-END READS") 

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

  #Additional usearch options related to alignment length
  temp_options <- paste0(' -fastq_minovlen ', dimsum_meta[["usearchMinovlen"]])
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options <- paste0(temp_options, " -xdrop_nw ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options <- paste0(temp_options, " -minhsp ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options <- paste0(temp_options, " -band ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 5 ){temp_options <- paste0(temp_options, " -hspw ", dimsum_meta[["usearchMinovlen"]])}

  #Run USEARCH on all fastq file pairs
  dimsum__status_message("Aligning paired-end FASTQ files with USEARCH:\n")
  dimsum__status_message(paste0(all_fastq, "\n"))
  dimsum__status_message("Processing...\n")
  #Trans library mode?
  if(dimsum_meta[["transLibrary"]]){
    for(i in 1:length(sample_names)){dimsum__status_message(paste0("\t", paste0(unlist(dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')]), collapse = "\t"), "\n"))}
    if(execute){
      # Setup cluster
      clust <- parallel::makeCluster(dimsum_meta[['numCores']])
      # make variables available to each core's workspace
      parallel::clusterExport(clust, list("dimsum_meta","usearch_outpath","sample_names","dimsum__concatenate_reads"), envir = environment())
      parallel::parSapply(clust,X = 1:length(sample_names), dimsum__usearch_trans_library_helper)
      parallel::stopCluster(clust)
    }
  #Single-end mode?
  }else if(!dimsum_meta[["paired"]]){
    for(i in 1:length(sample_names)){dimsum__status_message(paste0("\t", paste0(unique(unlist(dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')])), collapse = "\t"), "\n"))}
    if(execute){
      # Setup cluster
      clust <- parallel::makeCluster(dimsum_meta[['numCores']])
      # make variables available to each core's workspace
      parallel::clusterExport(clust, list("dimsum_meta","usearch_outpath","sample_names","dimsum__filter_single_end_reads"), envir = environment())
      parallel::parSapply(clust,X = 1:length(sample_names), dimsum_stage_usearch_single_end_library_helper)
      parallel::stopCluster(clust)
    }
  #Classic paired-end mode
  }else{
    for(i in 1:length(sample_names)){
      dimsum__status_message(paste0("\t", paste0(unlist(dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')]), collapse = "\t"), "\n"))
      #Check if this system command should be executed
      if(execute){
        temp_out <- system(paste0(
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
          dimsum_meta[['numCores']],
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
  #Delete files when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      temp_out <- mapply(file.remove, dimsum_meta_new[["deleteIntermediateFiles"]], MoreArgs = list(ignore.stdout = T, ignore.stderr = T))
      temp_out <- mapply(file.create, dimsum_meta_new[["touchIntermediateFiles"]], MoreArgs = list(ignore.stdout = T, ignore.stderr = T))
    }else{
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], 
        file.path(usearch_outpath, dir(usearch_outpath, "*.usearch$")))
    }
  }
  #Generate usearch report
  if(report){
    tryCatch({
      dimsum_meta_new_report <- dimsum__usearch_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
      }, error=function(e){
        dimsum__status_message("There were problems while running 'dimsum__usearch_report'\n")
        dimsum_meta_new_report <- dimsum_meta_new
        })
    return(dimsum_meta_new_report)
  }
  return(dimsum_meta_new)
}
