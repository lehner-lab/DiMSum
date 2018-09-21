
#dimsum_stage_filter
#
# Run fastx_collapser2AAtable.py on all input files.
#
# dimsum_meta: an experiment metadata object (required)
# filter_outpath: fastx_collapser2AAtable.py output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_filter <- function(
  dimsum_meta,
  filter_outpath,
  execute = TRUE
  ){
  #Create filter directory (if doesn't already exist)
  filter_outpath <- gsub("/$", "", filter_outpath)
  create_dimsum_dir(filter_outpath, execute = execute, message = "DiMSum STAGE 8: FILTER", overwrite_dir = FALSE)  
  #Construct filtered variant count table with corresponding amino acid sequences
  message("Constructing filtered variant count table with corresponding amino acid sequences")
  all_fasta <- file.path(dimsum_meta[["exp_design"]]$aligned_pair_unique_directory, dimsum_meta[['exp_design']]$aligned_pair_unique)
  print(all_fasta)
  message("Processing...")
  for(read_pair in dimsum_meta[['exp_design']]$aligned_pair_unique){
    print(read_pair)
    #Check if this system command should be executed
    if(execute){
      temp_out = system(paste0(
        "fastx_collapser2AAtable.py -i ", 
        file.path(dimsum_meta[["exp_design"]]$aligned_pair_unique_directory, read_pair), 
        " -o ", 
        file.path(filter_outpath, paste0(read_pair, ".tsv")), 
        " -n ", 
        nchar(dimsum_meta[['wildtypeSequence']]),
        " > ",
        file.path(filter_outpath, paste0(read_pair, '.tsv.stdout')),
        " 2> ",
        file.path(filter_outpath, paste0(read_pair, '.tsv.stderr'))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Filtered fasta filenames
  dimsum_meta_new[["exp_design"]]$aligned_pair_unique_tsv <- paste0(dimsum_meta_new[["exp_design"]]$aligned_pair_unique, ".tsv")
  dimsum_meta_new[['exp_design']]$aligned_pair_unique_tsv_directory <- filter_outpath
  #Get filter results for all samples
  filter_files <- file.path(dimsum_meta_new[['exp_design']]$aligned_pair_unique_tsv_directory, gsub('.tsv$', '.tsv.stdout', dimsum_meta_new[['exp_design']][,'aligned_pair_unique_tsv']))
  filter_list <- list()
  for(i in 1:length(filter_files)){
    temp_out <- system(paste0("cat ", filter_files[i]), intern=TRUE)
    filter_list[[i]] <- as.integer(rev(unlist(strsplit(temp_out[grep('Sequences failed', temp_out)], '\\t')))[1])
  }
  dimsum_meta_new[['exp_design']]$filter_incorrect_length = unlist(filter_list)
  return(dimsum_meta_new)
}
