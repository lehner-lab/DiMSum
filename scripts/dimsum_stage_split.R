
#dimsum_stage_split
#
# Split all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# split_outpath: split FASTQ output path (required)
# execute: whether or not to execute the system command (default: TRUE)
# save_workspace: whether or not to save the current experiment metadata object (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_split <- function(
  dimsum_meta,
  split_outpath,
  execute = TRUE,
  save_workspace = TRUE
  ){
  #Create/overwrite split directory (if executed)
  split_outpath <- gsub("/$", "", split_outpath)
  create_dimsum_dir(split_outpath, execute = execute, message = "DiMSum STAGE 4: SPLIT")  
  fastq_pair_list <- dimsum_meta[['exp_design']][,c('pair1', 'pair2')]
  rownames(fastq_pair_list) = 1:dim(fastq_pair_list)[1]
  #Split FASTQ files
  message("Splitting FASTQ files:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
  print(all_fastq)
  message("Processing...")
  for(pair_name in rownames(fastq_pair_list)){
    message(paste0("\t", fastq_pair_list[pair_name,]))
    #Check if this system command should be executed
    if(execute){
      temp_out = system(paste0(
        "fastq_splitter.py -i ", 
        file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][1]), 
        " -o ", 
        file.path(split_outpath, paste0(fastq_pair_list[pair_name,][1], ".split")), 
        " -c 3758096384",
        " > ",
        file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][1]), ".split.stdout")),
        " 2> ",
        file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][1]), ".split.stderr"))))
      num_records = as.integer(read.table(file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][1]), ".split.stdout"))))
      temp_out = system(paste0(
        "fastq_splitter.py -i ", 
        file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][2]), 
        " -o ", 
        file.path(split_outpath, paste0(fastq_pair_list[pair_name,][2], ".split")), 
        " -n ", 
        num_records,
        " >> ",
        file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][2]), ".split.stdout")),
        " 2>> ",
        file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][2]), ".split.stderr"))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  all_fastq <- file.path(dimsum_meta_new[["exp_design"]]$pair_directory, c(dimsum_meta_new[['exp_design']]$pair1))
  split_list <- list()
  for(f in all_fastq){
    num_splits <- length(list.files(split_outpath, pattern = basename(f)))
    split_list = append(split_list, num_splits)
  }
  dimsum_meta_new[["exp_design"]] = dimsum_meta_new[["exp_design"]][rep(1:length(all_fastq), times = unlist(split_list)),]
  temp_rownames = rownames(dimsum_meta_new[["exp_design"]])
  temp_suffix = rep('.split1', dim(dimsum_meta_new[["exp_design"]])[1])
  temp_suffix[grepl('\\.', temp_rownames)] = paste0('.split', as.integer(sapply(strsplit(temp_rownames[grepl('\\.', temp_rownames)], '\\.'), '[', 2))+1)
  dimsum_meta_new[["exp_design"]]$pair1 = paste0(dimsum_meta_new[["exp_design"]]$pair1, temp_suffix, '.fastq')
  dimsum_meta_new[["exp_design"]]$pair2 = paste0(dimsum_meta_new[["exp_design"]]$pair2, temp_suffix, '.fastq')
  dimsum_meta_new[["exp_design"]]$split = as.integer(gsub(".split", "", temp_suffix))
  dimsum_meta_new[['exp_design']]$pair_directory <- split_outpath
  #Save workspace
  if(save_workspace){save_metadata(dimsum_meta_new)}
  return(dimsum_meta_new)
}

