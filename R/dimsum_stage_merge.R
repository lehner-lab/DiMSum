
#' dimsum_stage_merge
#'
#' Merge all variant files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param merge_outpath merged variant data output path (required)
#' @param execute whether or not to execute the system command (default: TRUE)
#' @param report whether or not to generate final summary plots (default: TRUE)
#' @param report_outpath final summary report output path
#' @param save_workspace whether or not to save the current experiment metadata object (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
#' @import data.table
dimsum_stage_merge <- function(
  dimsum_meta,
  merge_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){
  #Create merge directory (if doesn't already exist)
  merge_outpath <- gsub("/$", "", merge_outpath)
  create_dimsum_dir(merge_outpath, execute = execute, message = "DiMSum STAGE 6: MERGE SAMPLE STATISTICS", overwrite_dir = FALSE)  
  #WT nucleotide sequence
  wt_NTseq <- tolower(dimsum_meta[['wildtypeSequence']])
  #WT AA sequence
  wt_AAseq <- paste0(seqinr::translate(strsplit(wt_NTseq,split="")[[1]]),collapse="")
  #Initialise merged variant data object
  variant_data <- NULL
  #Nucleotide mutation count dictionaries
  nuc_subst_dict <- list()
  nuc_indel_dict <- list()
  #AA mutation count dictionaries
  aa_subst_dict <- list()
  aa_indel_dict <- list()
  
  #Load variant count files
  message("Loading variant count files:")
  all_count <- file.path(dimsum_meta[['exp_design']][,"aligned_pair_unique_directory"], dimsum_meta[['exp_design']][,"aligned_pair_unique"])
  print(all_count)
  message("Processing...")
  for(count_file in all_count){
    message(paste0("\t", count_file))
    #Check if this code should be executed
    if(execute){
      file_id <- gsub('.usearch.unique', '', basename(count_file))
      #Initialise count table
      count_dt <- data.table(
        nt_seq = character(),
        count = numeric(),
        aa_seq = character(),
        Nins_aa = integer(),
        Ndel_aa = integer(),
        Nsub_aa = integer(),
        Nmut_aa = integer(),
        Nins_nt = integer(),
        Ndel_nt = integer(),
        Nsub_nt = integer(),
        Nmut_nt = integer())
      #Load fasta file (if exists)
      if(!file.exists(count_file)){
        warning("File does not exist.", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
        #Save nucleotide mutation distribution (with/without indels)
        nuc_subst_dict[[count_file]] <- NA
        nuc_indel_dict[[count_file]] <- NA
        #Save amino acid mutation distribution (with/without indels)
        aa_subst_dict[[count_file]] <- NA
        aa_indel_dict[[count_file]] <- NA
      }else{
        rfa <- ShortRead::readFasta(count_file)
        #Create variant count table with nucleotide and amino acid sequence
        count_dt <- data.table(nt_seq = tolower(as.character(ShortRead::sread(rfa))))
        count_dt[, count := as.numeric(sapply(strsplit(as.character(ShortRead::id(rfa)), "-"), "[", 2))]
        suppressWarnings(count_dt[, aa_seq := as.character(Biostrings::translate(ShortRead::sread(rfa)))])
        #Calculate number of aa mutations (insertions, deletions, substitutions)
        mut_counts <- attr(utils::adist(count_dt[,aa_seq], wt_AAseq, counts = T), "counts")
        count_dt[,Nins_aa := mut_counts[,1,2]]
        count_dt[,Ndel_aa := mut_counts[,1,1]]
        count_dt[,Nsub_aa := mut_counts[,1,3]]
        count_dt[,Nmut_aa := Nins_aa+Ndel_aa+Nsub_aa]
        #Calculate number of nucleotide mutations (insertions, deletions, substitutions)
        mut_counts <- attr(utils::adist(count_dt[,nt_seq], wt_NTseq, counts = T), "counts")
        count_dt[,Nins_nt := mut_counts[,1,2]]
        count_dt[,Ndel_nt := mut_counts[,1,1]]
        count_dt[,Nsub_nt := mut_counts[,1,3]]
        count_dt[,Nmut_nt := Nins_nt+Ndel_nt+Nsub_nt]
        #Save nucleotide mutation distribution (with/without indels)
        nuc_subst_dict[[count_file]] <- tapply(count_dt[Nsub_nt==Nmut_nt,count], count_dt[Nsub_nt==Nmut_nt,Nsub_nt], sum)
        nuc_indel_dict[[count_file]] <- tapply(count_dt[Nsub_nt!=Nmut_nt,count], count_dt[Nsub_nt!=Nmut_nt, .(Nindels_nt = Nins_nt + Ndel_nt)][,Nindels_nt], sum)
        #Save amino acid mutation distribution (with/without indels)
        aa_subst_dict[[count_file]] <- tapply(count_dt[Nsub_aa==Nmut_aa,count], count_dt[Nsub_aa==Nmut_aa,Nsub_aa], sum)
        aa_indel_dict[[count_file]] <- tapply(count_dt[Nsub_aa!=Nmut_aa,count], count_dt[Nsub_aa!=Nmut_aa, .(Nindels_aa = Nins_aa + Ndel_aa)][,Nindels_aa], sum)
      }
      #First file loaded
      if(is.null(variant_data)){
        variant_data <- count_dt
      #Not first file loaded (merge with previous data)
      }else{
        variant_data = merge(variant_data,count_dt,by=names(variant_data)[-grep("count", names(variant_data))],all = T)
      }
      names(variant_data)[names(variant_data)=="count"] = paste0(file_id, "_count")
    }
  }
  #Check if this code should be executed
  if(execute){
    #Save mutation statistics dictionaries
    mutation_stats_dicts <- list(
      "nuc_subst_dict" = nuc_subst_dict, 
      "nuc_indel_dict" = nuc_indel_dict, 
      "aa_subst_dict" = aa_subst_dict,
      "aa_indel_dict" = aa_indel_dict)
    save(mutation_stats_dicts, file = file.path(dimsum_meta[["tmp_path"]], "mutation_stats_dicts.RData"))
    #Replace NA counts with zeros
    variant_data[is.na(variant_data)] <- 0
    #Indicate WT sequence
    variant_data[nt_seq == wt_NTseq,WT := TRUE]
    #Indicate STOPs
    variant_data[,STOP := ifelse(length(grep(aa_seq,pattern="\\*"))==1,TRUE,FALSE),aa_seq]
    #Merge split counts
    split_base <- unique(sapply(strsplit(colnames(variant_data)[grep('_count', colnames(variant_data))], "_split"), '[', 1))
    variant_data_merge <- sum_datatable_columns(dt=variant_data, column_patterns=split_base, suffix="_count")
    #Merge technical counts
    split_base <- unique(sapply(strsplit(colnames(variant_data_merge)[grep('_count', colnames(variant_data_merge))], "_t"), '[', 1))
    variant_data_merge <- sum_datatable_columns(dt=variant_data_merge, column_patterns=split_base, suffix="_count")
    #Save merged variant data
    message("Saving merged variant data...")
    save(variant_data_merge, file=file.path(merge_outpath, paste0(dimsum_meta[["project_name"]], '_variant_data_merge.RData')))
    message("Done")
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Merged variant data path
  dimsum_meta_new[["variant_data_merge_path"]] <- file.path(merge_outpath, paste0(dimsum_meta[["project_name"]], '_variant_data_merge.RData')) 
  #Load mutation statistics
  load(file.path(dimsum_meta[["tmp_path"]], "mutation_stats_dicts.RData"))
  #AA mutation results
  dimsum_meta_new[["aa_subst_counts"]] <- mutation_stats_dicts[["aa_subst_dict"]]
  dimsum_meta_new[["aa_indel_counts"]] <- mutation_stats_dicts[["aa_indel_dict"]]
  #Nucleotide mutation results
  dimsum_meta_new[["nuc_subst_counts"]] <- mutation_stats_dicts[["nuc_subst_dict"]]
  dimsum_meta_new[["nuc_indel_counts"]] <- mutation_stats_dicts[["nuc_indel_dict"]]
  #Generate merge report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_merge_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
    dimsum_stage_diagnostics_report(variant_data = file.path(merge_outpath, paste0(dimsum_meta[["project_name"]], '_variant_data_merge.RData')), report_outpath = report_outpath)
    #Save workspace
    if(save_workspace){save_metadata(dimsum_meta_new_report)}
    return(dimsum_meta_new_report)
  }
  #Save workspace
  if(save_workspace){save_metadata(dimsum_meta_new)}
  return(dimsum_meta_new)
}
