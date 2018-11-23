
#dimsum_stage_merge_indels
#
# Merge all indel variant files.
#
# dimsum_meta: an experiment metadata object (required)
# merge_outpath: merged variant data output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_merge_indels <- function(
  dimsum_meta,
  merge_outpath,
  execute = TRUE
  ){
  #Create merge directory (if doesn't already exist)
  merge_outpath <- gsub("/$", "", merge_outpath)
  create_dimsum_dir(merge_outpath, execute = execute, message = "DiMSum STAGE 10: MERGE INDELS", overwrite_dir = FALSE)  
  #WT nucleotide sequence
  wt_NTseq <- tolower(dimsum_meta[['wildtypeSequence']])
  #WT AA sequence
  wt_AAseq <- paste0(seqinr::translate(strsplit(wt_NTseq,split="")[[1]]),collapse="")
  #Initialise dictionary of AA mappings
  aa_dict <- list()
  #Initialise merged variant data object
  variant_data <- NULL
  
  #Load filtered variant count files
  message("Loading filtered indel variant count files:")
  all_count <- file.path(dimsum_meta[['exp_design']]$aligned_pair_unique_indel_tsv_directory, dimsum_meta[['exp_design']]$aligned_pair_unique_indel_tsv)
  print(all_count)
  message("Processing...")
  for(count_file in all_count){
    print(count_file)
    #Check if this code should be executed
    if(execute){
      file_id <- gsub('.usearch.unique.indels.tsv', '', basename(count_file))
      #Load file
      temp_file <- fread(count_file, header = T, sep="\t", stringsAsFactors = F)
      #Calculate number of aa mutations
      temp_file[,Nmut_aa := adist(aa_seq,wt_AAseq)]
      #Calculate number of nucleotide mutations
      temp_file[,Nmut_nt := adist(nt_seq,wt_NTseq)]
      #Save AA mapping to dictionary
      temp_dict <- as.list(temp_file$aa_seq)
      names(temp_dict) <- temp_file$nt_seq
      aa_dict <- c(aa_dict, temp_dict)
      aa_dict <- aa_dict[unique(names(aa_dict))]
      #First file loaded
      if(count_file == all_count[1]){
        variant_data <- temp_file
        names(variant_data)[grep(names(variant_data),pattern="count")] = paste0(file_id, "_count")
      #Not first file loaded (merge with previous data)
      }else{
        variant_data = merge(variant_data,temp_file[,.(nt_seq,.SD),,.SDcols=names(temp_file)[grep(names(temp_file),pattern="count")]],by="nt_seq",all = T)
        names(variant_data)[grep(names(variant_data),pattern=".SD")] = paste0(file_id, "_count")
      }
    }
  }
  #Check if this code should be executed
  if(execute){
    #Add AA sequence to merged variant data
    variant_data$aa_seq <- unlist(aa_dict[variant_data$nt_seq])
    #Calculate number of AA mutations of merged variant data
    variant_data[,Nmut_aa := adist(aa_seq,wt_AAseq)]
    #Replace NA counts with zeros
    variant_data[is.na(variant_data)] <- 0
    #Indicate WT sequence
    variant_data[nt_seq == wt_NTseq,WT := TRUE]
    #Indicate STOPs
    variant_data[,STOP := ifelse(length(grep(aa_seq,pattern="\\*"))==1,TRUE,FALSE),aa_seq]
    #Calculate number of nucleotide mutations of merged variant data
    variant_data[,Nmut_nt := adist(nt_seq,wt_NTseq)]
    #Merge split counts
    split_base <- unique(sapply(strsplit(colnames(variant_data)[grep('_count', colnames(variant_data))], "_split"), '[', 1))
    variant_data_merge_indels <- sum_datatable_columns(dt=variant_data, column_patterns=split_base, suffix="_count")
    #Merge technical counts
    split_base <- unique(sapply(strsplit(colnames(variant_data_merge_indels)[grep('_count', colnames(variant_data_merge_indels))], "_t"), '[', 1))
    variant_data_merge_indels <- sum_datatable_columns(dt=variant_data_merge_indels, column_patterns=split_base, suffix="_count")
    #Save merged variant data
    message("Saving merged variant data...")
    save(variant_data_merge_indels, file=file.path(merge_outpath, paste0(dimsum_meta[["project_name"]], '_variant_data_merge_indels.RData')))
    message("Done")
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Merged variant data path
  dimsum_meta_new[["variant_data_merge_indels_path"]] <- file.path(merge_outpath, paste0(dimsum_meta[["project_name"]], '_variant_data_merge_indels.RData')) 
  return(dimsum_meta_new)
}
