
#' dimsum_stage_split
#'
#' Split all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param split_outpath split FASTQ output path (required)
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_split <- function(
  dimsum_meta,
  split_outpath,
  save_workspace = TRUE
  ){
  #Whether or not to execute the system command
  this_stage <- 2
  execute <- (dimsum_meta[["startStage"]] <= this_stage & (dimsum_meta[["stopStage"]] == 0 | dimsum_meta[["stopStage"]] >= this_stage))
  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}
  #Create/overwrite split directory (if executed)
  split_outpath <- gsub("/$", "", split_outpath)
  dimsum__create_dir(split_outpath, execute = execute, message = "SPLIT FASTQ FILES")  
  fastq_pair_list <- unique(dimsum_meta[['exp_design']][,c('pair1', 'pair2')])
  rownames(fastq_pair_list) <- 1:dim(fastq_pair_list)[1]
  #Split FASTQ files
  message("Splitting FASTQ files:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], unique(c(dimsum_meta[['exp_design']][,"pair1"], dimsum_meta[['exp_design']][,"pair2"])))
  print(all_fastq)
  message("Processing...")
  message(paste0("\t", basename(all_fastq), "\n"))
  #Check if this system command should be executed
  if(execute){
    dimsum_stage_split_helper <- function(
      i
      ){
      num_records <- dimsum__fastq_splitter(
        inputFile = file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[i,][1]),
        outputFilePrefix = file.path(split_outpath, paste0(fastq_pair_list[i,][1], ".split")),
        chunkSize = dimsum_meta[["splitChunkSize"]])
      if(dimsum_meta[["paired"]]){
        num_records <- dimsum__fastq_splitter(
          inputFile = file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[i,][2]),
          outputFilePrefix = file.path(split_outpath, paste0(fastq_pair_list[i,][2], ".split")),
          numRecords = num_records)
      }
    }
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","fastq_pair_list","split_outpath","dimsum__fastq_splitter","dimsum__writeFastq"), envir = environment())
    parallel::parSapply(clust,X = 1:nrow(fastq_pair_list), dimsum_stage_split_helper)
    parallel::stopCluster(clust)
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  all_fastq <- file.path(dimsum_meta_new[["exp_design"]][,"pair_directory"], c(dimsum_meta_new[['exp_design']][,"pair1"]))
  split_list <- list()
  for(f in all_fastq){
    num_splits <- length(list.files(split_outpath, pattern = basename(f)))
    split_list <- append(split_list, num_splits)
  }
  dimsum_meta_new[["exp_design"]] <- dimsum_meta_new[["exp_design"]][rep(1:length(all_fastq), times = unlist(split_list)),]
  temp_rownames <- rownames(dimsum_meta_new[["exp_design"]])
  temp_suffix <- rep('.split1', dim(dimsum_meta_new[["exp_design"]])[1])
  temp_suffix[grepl('\\.', temp_rownames)] <- paste0('.split', as.integer(sapply(strsplit(temp_rownames[grepl('\\.', temp_rownames)], '\\.'), '[', 2))+1)
  dimsum_meta_new[["exp_design"]][,"pair1"] <- paste0(dimsum_meta_new[["exp_design"]][,"pair1"], temp_suffix, '.fastq')
  dimsum_meta_new[["exp_design"]][,"pair2"] <- paste0(dimsum_meta_new[["exp_design"]][,"pair2"], temp_suffix, '.fastq')
  dimsum_meta_new[["exp_design"]][,"split"] <- as.integer(gsub(".split", "", temp_suffix))
  dimsum_meta_new[['exp_design']][,"pair_directory"] <- split_outpath
  #Delete file contents when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      temp_out <- mapply(system, dimsum_meta_new[["deleteIntermediateFiles"]], MoreArgs = list(ignore.stdout = T, ignore.stderr = T))
    }else{
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], 
        paste0("> ", file.path(dimsum_meta_new[['exp_design']][,"pair_directory"], c(dimsum_meta_new[['exp_design']][,"pair1"], dimsum_meta_new[['exp_design']][,"pair2"]))))
    }
  }
  return(dimsum_meta_new)
}

