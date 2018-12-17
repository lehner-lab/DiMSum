
#dimsum_stage_demultiplex
#
# Run demultiplexing on all fastq files using cutadapt.
#
# dimsum_meta: an experiment metadata object (required)
# demultiplex_outpath: demultiplex output path (required)
# execute: whether or not to execute the system command (default: TRUE)
# save_workspace: whether or not to save the current experiment metadata object (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_demultiplex <- function(
  dimsum_meta,
  demultiplex_outpath,
  execute = TRUE,
  save_workspace = TRUE
  ){
  #Create/overwrite demultiplex directory (if executed)
  demultiplex_outpath <- gsub("/$", "", demultiplex_outpath)
  create_dimsum_dir(demultiplex_outpath, execute = execute, message = "DiMSum STAGE 1: DEMULTIPLEX")  
  #Demultiplex parameters specified?
  if( !'barcode_design' %in% names(dimsum_meta) ){
    message("Skipping this stage (assuming all fastq files already demultiplexed)")
    #Save workspace
    if(save_workspace){save_metadata(dimsum_meta)}
    return(dimsum_meta)
  }
  fastq_pair_list <- unique(dimsum_meta[['barcode_design']][,c('pair1', 'pair2')])
  rownames(fastq_pair_list) = 1:dim(fastq_pair_list)[1]
  #Reformat barcode files for cutadapt (and copy and rename FASTQ files if necessary)
  for(pair_name in rownames(fastq_pair_list)){
    #Check if this system command should be executed
    if(execute){
      #
      temp_design <- dimsum_meta[['barcode_design']][dimsum_meta[['barcode_design']][,'pair1']==fastq_pair_list[pair_name,'pair1'],]
      write(
        x = c(rbind(paste0('>', temp_design[,"new_pair_prefix"]), paste0('^', temp_design[,"barcode"]))), 
        file = file.path(demultiplex_outpath, paste0('demultiplex_barcode-file_', pair_name, '.fasta')), 
        sep="\n")
      #Check if file extension incompatible with cutadapt (i.e. NOT ".fastq" or ".fastq.gz")
      if(dimsum_meta[["fastq_file_extension"]]!=".fastq"){
        #Copy FASTQ files to temp directory and format extension
        new_fastq_name1 <- gsub(paste0(dimsum_meta[["fastq_file_extension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), c(".fastq", ".fastq.gz")[as.numeric(dimsum_meta[["gzipped"]])+1], fastq_pair_list[pair_name,][1])
        new_fastq_name2 <- gsub(paste0(dimsum_meta[["fastq_file_extension"]], c("$", ".gz$")[as.numeric(dimsum_meta[["gzipped"]])+1]), c(".fastq", ".fastq.gz")[as.numeric(dimsum_meta[["gzipped"]])+1], fastq_pair_list[pair_name,][2])
        temp_out = system(paste0(
          "cp ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][1]),
          " ",
          file.path(demultiplex_outpath, new_fastq_name1)))
        temp_out = system(paste0(
          "cp ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"][1], fastq_pair_list[pair_name,][2]),
          " ",
          file.path(demultiplex_outpath, new_fastq_name2)))
        #Update names in list
        fastq_pair_list[pair_name,][1] <- new_fastq_name1
        fastq_pair_list[pair_name,][2] <- new_fastq_name2
      }
    }
  }
  #Check if file extension incompatible with cutadapt (i.e. NOT ".fastq" or ".fastq.gz")
  if(dimsum_meta[["fastq_file_extension"]]!=".fastq"){
    #Update FASTQ file directory in list
    dimsum_meta[["exp_design"]][,"pair_directory"] <- demultiplex_outpath
  }
  #Demultiplex FASTQ files
  message("Demultiplexing FASTQ files with cutadapt:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]][,"pair_directory"], c(dimsum_meta[['barcode_design']][,"pair1"], dimsum_meta[['barcode_design']][,"pair2"]))
  print(unique(all_fastq))
  message("Processing...")
  for(pair_name in rownames(fastq_pair_list)){
    #TODO: saber binary path specifiable on commandline?
    message(paste0("\t", fastq_pair_list[pair_name,]))
    #Check if this system command should be executed
    if(execute){
      #Demultiplex using cutadapt
      temp_out = system(paste0(
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
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']][,"pair_directory"] <- demultiplex_outpath
  dimsum_meta_new[["fastq_file_extension"]] <- ".fastq"
  dimsum_meta_new[["gzipped"]] <- FALSE
  #Save workspace
  if(save_workspace){save_metadata(dimsum_meta_new)}
  return(dimsum_meta_new)
}

