
#' dimsum_stage_merge
#'
#' Merge all variant files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param merge_outpath merged variant data output path (required)
#' @param report whether or not to generate final summary plots (default: TRUE)
#' @param report_outpath final summary report output path
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
#' @import data.table
dimsum_stage_merge <- function(
  dimsum_meta,
  merge_outpath,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 4
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Create merge directory (if doesn't already exist)
  merge_outpath <- gsub("/$", "", merge_outpath)
  dimsum__create_dir(merge_outpath, execute = execute, message = "DiMSum STAGE 4: PROCESS VARIANT SEQUENCES", overwrite_dir = FALSE)  

  #Input files
  if(!is.null(dimsum_meta[["countPath"]])){
    all_count <- dimsum_meta[["countPath"]]
  }else{
    all_count <- unique(file.path(dimsum_meta[['exp_design']][,"aligned_pair_unique_directory"], dimsum_meta[['exp_design']][,"aligned_pair_unique"]))
  }

  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = all_count,
    stage_number = this_stage,
    execute = execute)

  #WT nucleotide sequence
  wt_ntseq <- tolower(dimsum_meta[['wildtypeSequence']])
  #WT AA sequence
  wt_AAseq <- paste0(seqinr::translate(strsplit(wt_ntseq,split="")[[1]]),collapse="")
  #Initialise merged variant data object
  variant_data <- NULL
  
  #Load variant count files
  dimsum__status_message("Loading variant count files:\n")
  dimsum__status_message(paste0(all_count, "\n"))
  dimsum__status_message("Processing...\n")
  for(count_file in all_count){
    dimsum__status_message(paste0("\t", basename(count_file), "\n"))
    #Check if this code should be executed
    if(execute){
      
      #load variant data from user-specified count file
      if(!is.null(dimsum_meta[["countPath"]])){
        variant_data <- data.table::fread(count_file)
        variant_data <- dimsum__check_countfile(dimsum_meta = dimsum_meta, input_dt = variant_data)
        break
      }

      #load variant data from dimsum variant count files
      file_id <- gsub('.usearch.unique', '', basename(count_file))
      #Initialise count table
      count_dt <- data.table(
        nt_seq = character(),
        count = numeric())
      #Load fasta file (if exists)
      if(!file.exists(count_file)){
        warning("dimsum_stage_merge.R: File does not exist.", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
      }else{
        #Create variant count table with nucleotide and amino acid sequence
        temp_dt <- fread(count_file)
        if(dim(temp_dt)[1]!=0){count_dt <- temp_dt}
        names(count_dt) <- c("nt_seq", "count")
        count_dt[, nt_seq := tolower(nt_seq)]
      }
      #First file loaded
      if(is.null(variant_data)){
        variant_data <- count_dt
      #Not first file loaded (merge with previous data)
      }else{
        variant_data <- merge(variant_data,count_dt,by=names(variant_data)[-grep("count", names(variant_data))],all = T)
      }
      names(variant_data)[names(variant_data)=="count"] = paste0(file_id, "_count")
    }
  }
  #Check if this code should be executed
  if(execute){
    #Replace NA counts with zeros
    variant_data[is.na(variant_data)] <- 0

    #Merge counts from technical replicates
    split_base <- unique(sapply(strsplit(colnames(variant_data)[grep('_count', colnames(variant_data))], "_t"), '[', 1))
    variant_data_merge <- dimsum__sum_datatable_columns(dt=variant_data, column_patterns=split_base, suffix="_count")

    #Process variants to separate into indel, rejected and retained variants data.tables
    variant_list <- dimsum__process_merged_variants(dimsum_meta = dimsum_meta, input_dt = variant_data_merge)

    #Save merged variant data
    dimsum__status_message("Saving merged variant data...\n")
    write.table(variant_list[["nobarcode_variants"]], file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_nobarcode_variant_data_merge.tsv')), sep = "\t", quote = F, row.names = F)
    write.table(variant_list[["indel_variants"]], file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_indel_variant_data_merge.tsv')), sep = "\t", quote = F, row.names = F)
    write.table(variant_list[["rejected_variants"]], file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_rejected_variant_data_merge.tsv')), sep = "\t", quote = F, row.names = F)
    variant_data_merge <- variant_list[["retained_variants"]]
    save(variant_data_merge, 
      file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_variant_data_merge.RData')),
      version = 2)
    write.table(variant_list[["retained_variants"]], file=file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_variant_data_merge.tsv')), sep = "\t", quote = F, row.names = F)
    dimsum__status_message("Done\n")
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Merged variant data path
  dimsum_meta_new[["variant_data_merge_path"]] <- file.path(merge_outpath, paste0(dimsum_meta[["projectName"]], '_variant_data_merge.RData')) 
  #Load previously saved mutation statistics
  load(file.path(dimsum_meta[["tmp_path"]], "mutation_stats_dicts.RData"))
  #Add to metadata
  dimsum_meta_new[["aa_subst_counts"]] <- mutation_stats_dicts[["aa_subst_dict"]]
  dimsum_meta_new[["nuc_subst_counts"]] <- mutation_stats_dicts[["nuc_subst_dict"]]
  dimsum_meta_new[["nuc_mxsub_counts"]] <- mutation_stats_dicts[["nuc_mxsub_dict"]]
  dimsum_meta_new[["nuc_tmsub_counts"]] <- mutation_stats_dicts[["nuc_tmsub_dict"]]
  dimsum_meta_new[["nuc_frbdn_counts"]] <- mutation_stats_dicts[["nuc_frbdn_dict"]]
  dimsum_meta_new[["nuc_const_counts"]] <- mutation_stats_dicts[["nuc_const_dict"]]
  dimsum_meta_new[["nuc_indel_counts"]] <- mutation_stats_dicts[["nuc_indel_dict"]]
  dimsum_meta_new[["nuc_nbarc_counts"]] <- mutation_stats_dicts[["nuc_nbarc_dict"]]
  #Delete files when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      temp_out <- mapply(file.remove, dimsum_meta_new[["deleteIntermediateFiles"]], MoreArgs = list(ignore.stdout = T, ignore.stderr = T))
      temp_out <- mapply(file.create, dimsum_meta_new[["touchIntermediateFiles"]], MoreArgs = list(ignore.stdout = T, ignore.stderr = T))
    }
  }
  #Generate merge report
  if(report){
    tryCatch({
      dimsum_meta_new_report <- dimsum__merge_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
      }, error=function(e){
        dimsum__status_message("There were problems while running 'dimsum__merge_report'\n")
        dimsum_meta_new_report <- dimsum_meta_new
        })
    tryCatch({
      dimsum__diagnostics_report(dimsum_meta = dimsum_meta_new_report, report_outpath = report_outpath)
      }, error=function(e){
        dimsum__status_message("There were problems while running 'dimsum__diagnostics_report'\n")
        })
    return(dimsum_meta_new_report)
  }
  return(dimsum_meta_new)
}
