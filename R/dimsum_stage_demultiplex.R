
#' dimsum_stage_demultiplex
#'
#' Run demultiplexing on all fastq files using cutadapt.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param demultiplex_outpath demultiplex output path (required)
#' @param execute whether or not to execute the system command (default: TRUE)
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_demultiplex <- function(
  dimsum_meta,
  demultiplex_outpath,
  execute = TRUE,
  save_workspace = TRUE
  ){
  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}
  #Create/overwrite demultiplex directory (if executed)
  demultiplex_outpath <- gsub("/$", "", demultiplex_outpath)
  dimsum__create_dir(demultiplex_outpath, execute = execute, message = "DiMSum STAGE 1: DEMULTIPLEX READS")  
  #Demultiplex parameters specified?
  if( !'barcode_design' %in% names(dimsum_meta) ){
    message("Skipping this stage (assuming all fastq files already demultiplexed)")
    return(dimsum_meta)
  }
  fastq_pair_list <- unique(dimsum_meta[['barcode_design']][,c('pair1', 'pair2')])
  rownames(fastq_pair_list) = 1:dim(fastq_pair_list)[1]
  #Reformat barcode files for cutadapt (and copy and rename FASTQ files if necessary)
  #Check if this system command should be executed
  if(execute){
    dimsum_stage_demultiplex_cp_helper <- function(
      i
      ){
      pair_name <- rownames(fastq_pair_list)[i]
      temp_design <- dimsum_meta[['barcode_design']][dimsum_meta[['barcode_design']][,'pair1']==fastq_pair_list[pair_name,'pair1'],]
      write(
        x = c(rbind(paste0('>', temp_design[,"new_pair_prefix"]), paste0('^', temp_design[,"barcode"]))), 
        file = file.path(demultiplex_outpath, paste0('demultiplex_barcode-file_', pair_name, '.fasta')), 
        sep="\n")
      #Check if file extension incompatible with cutadapt (i.e. NOT ".fastq" or ".fastq.gz")
      if(dimsum_meta[["fastqFileExtension"]]!=".fastq"){
        #Copy FASTQ files to temp directory and format extension
        new_fastq_name1 <- gsub(paste0(dimsum_meta[["fastqFileExtension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), c(".fastq", ".fastq.gz")[as.numeric(dimsum_meta[["gzipped"]])+1], fastq_pair_list[pair_name,][1])
        new_fastq_name2 <- gsub(paste0(dimsum_meta[["fastqFileExtension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), c(".fastq", ".fastq.gz")[as.numeric(dimsum_meta[["gzipped"]])+1], fastq_pair_list[pair_name,][2])
        temp_out <- system(paste0(
          "cp ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][1]),
          " ",
          file.path(demultiplex_outpath, new_fastq_name1)))
        #If second read in pair exists
        if(dimsum_meta[["paired"]]){
          temp_out <- system(paste0(
            "cp ",
            file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][2]),
            " ",
            file.path(demultiplex_outpath, new_fastq_name2)))
        }
      }
    }
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","fastq_pair_list","demultiplex_outpath"), envir = environment())
    parallel::parSapply(clust,X = 1:nrow(fastq_pair_list), dimsum_stage_demultiplex_cp_helper)
    parallel::stopCluster(clust)
  }
  #Update names in list
  for(pair_name in rownames(fastq_pair_list)){
    #Check if this system command should be executed
    if(execute){
      #Check if file extension incompatible with cutadapt (i.e. NOT ".fastq" or ".fastq.gz")
      if(dimsum_meta[["fastqFileExtension"]]!=".fastq"){
        #New FASTQ file names
        new_fastq_name1 <- gsub(paste0(dimsum_meta[["fastqFileExtension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), c(".fastq", ".fastq.gz")[as.numeric(dimsum_meta[["gzipped"]])+1], fastq_pair_list[pair_name,][1])
        new_fastq_name2 <- gsub(paste0(dimsum_meta[["fastqFileExtension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), c(".fastq", ".fastq.gz")[as.numeric(dimsum_meta[["gzipped"]])+1], fastq_pair_list[pair_name,][2])
        #Update names in list
        fastq_pair_list[pair_name,][1] <- new_fastq_name1
        fastq_pair_list[pair_name,][2] <- new_fastq_name2
      }
    }
  }
  #Check if file extension incompatible with cutadapt (i.e. NOT ".fastq" or ".fastq.gz")
  if(dimsum_meta[["fastqFileExtension"]]!=".fastq"){
    #Update FASTQ file directory in list
    dimsum_meta[["exp_design"]][,"pair_directory"] <- demultiplex_outpath
  }
  #Demultiplex FASTQ files
  message("Demultiplexing FASTQ files with cutadapt:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]][,"pair_directory"], unique(c(dimsum_meta[['barcode_design']][,"pair1"], dimsum_meta[['barcode_design']][,"pair2"])))
  print(unique(all_fastq))
  message("Processing...")
  for(i in 1:dim(fastq_pair_list)[1]){message(paste0("\t", unique(unlist(fastq_pair_list[i,]))))}
  #Check if this system command should be executed
  if(execute){
    dimsum_stage_demultiplex_helper <- function(
      i
      ){
      pair_name <- rownames(fastq_pair_list)[i]
      #Demultiplex using cutadapt
      if(dimsum_meta[["paired"]]){
        temp_out <- system(paste0(
          "cutadapt",
          " -g file:",
          file.path(demultiplex_outpath, paste0('demultiplex_barcode-file_', pair_name, '.fasta')),
          " -G file:",
          file.path(demultiplex_outpath, paste0('demultiplex_barcode-file_', pair_name, '.fasta')),
          " -e ",
          as.character(dimsum_meta[["barcodeErrorRate"]]),
          " --no-indels ",
          " --untrimmed-output ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.unknown.fastq")),
          " --untrimmed-paired-output ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][2], ".demultiplex.unknown.fastq")),
          " -o ",
          file.path(demultiplex_outpath, "{name}1.fastq"),
          " -p ",
          file.path(demultiplex_outpath, "{name}2.fastq"),
          " ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][1]),
          " ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][2]),
          " > ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.stdout")),
          " 2> ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.stderr"))))
      }else{
        temp_out <- system(paste0(
          "cutadapt",
          " -g file:",
          file.path(demultiplex_outpath, paste0('demultiplex_barcode-file_', pair_name, '.fasta')),
          " -e ",
          as.character(dimsum_meta[["barcodeErrorRate"]]),
          " --no-indels ",
          " --untrimmed-output ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.unknown.fastq")),
          " -o ",
          file.path(demultiplex_outpath, "{name}1.fastq"),
          " ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][1]),
          " > ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.stdout")),
          " 2> ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.stderr"))))        
      }
    }
    # Setup cluster
    clust <- parallel::makeCluster(dimsum_meta[['numCores']])
    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("dimsum_meta","demultiplex_outpath","fastq_pair_list"), envir = environment())
    parallel::parSapply(clust,X = 1:nrow(fastq_pair_list), dimsum_stage_demultiplex_helper)
    parallel::stopCluster(clust)
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']][,"pair_directory"] <- demultiplex_outpath
  dimsum_meta_new[["fastqFileExtension"]] <- ".fastq"
  dimsum_meta_new[["gzipped"]] <- FALSE
  return(dimsum_meta_new)
}

