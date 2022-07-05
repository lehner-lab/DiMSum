
#' dimsum__process_merged_variants
#'
#' Process nucleotide variants to remove indels, internal constant region variants, not permitted sequences and variants with too many substitutions.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt output path for plots and saved objects (required)
#'
#' @return A list of data.tables containing indel variants, rejected variants (internal constant region variants, not permitted sequences and variants with too many substititions) and retained variants
#' @export
#' @import data.table
dimsum__process_merged_variants <- function(
  dimsum_meta,
  input_dt
  ){

  dimsum__status_message("Processing merged variants...\n")

  variant_dt <- input_dt

  #WT nucleotide sequence
  wt_ntseq <- tolower(dimsum_meta[['wildtypeSequence']])
  #WT AA sequence
  wt_AAseq <- paste0(seqinr::translate(strsplit(wt_ntseq,split="")[[1]]),collapse="")

  #Nucleotide mutation count dictionaries
  nuc_nbarc_dict <- list()
  nuc_indel_dict <- list()
  nuc_const_dict <- list()
  nuc_frbdn_dict <- list()
  nuc_tmsub_dict <- list()
  nuc_mxsub_dict <- list()
  nuc_indrt_dict <- list()
  nuc_subst_dict <- list()
  aa_subst_dict <- list()

  ### Reverse complement (if necessary)
  ###########################

  if(dimsum_meta[["reverseComplement"]]){
    variant_dt[,nt_seq := tolower(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nt_seq))))]
  }

  ### Debarcode variants (if barcode identity file supplied)
  ###########################

  variant_dt[, barcode_valid := TRUE]

  if(!is.null(dimsum_meta[["barcodeIdentityPath"]])){
    variant_dt <- dimsum__debarcode_variants(
      input_dt = variant_dt,
      barcode_path = dimsum_meta[["barcodeIdentityPath"]])
  }

  #No barcode variant count statistics
  nuc_nbarc_dict <- as.list(apply(
    variant_dt[barcode_valid==F,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]], 
    2, sum))

  #No barcode variants
  nobarcode_variant_dt <- data.table::copy(variant_dt[barcode_valid==F,.SD,,.SDcols = c(
    "nt_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid")])

  ### Update remaining variants
  ###########################

  #Remaining variants with valid barcodes
  variant_dt <- variant_dt[barcode_valid==T,.SD,,.SDcols = c(
    "nt_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid")]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No valid variant barcodes found"), call. = FALSE)
  }

  #Define indel variants as those with length different from WT
  variant_dt[, indel := F]
  variant_dt[nchar(nt_seq)!=nchar(wt_ntseq), indel := T]

  #Check if substitution variants exist (if indels discarded)
  if(variant_dt[indel==F,.N]==0 & !dimsum_meta[["indels"]]){
    stop(paste0("Cannot proceed with variant processing: No substitution variants found"), call. = FALSE)
  }

  #Calculate nucleotide hamming distance
  variant_dt[indel==F, Nham_nt := mapply(dimsum__hamming_distance, nt_seq, wt_ntseq)]

  #Determine amino acid sequence
  suppressWarnings(variant_dt[, aa_seq := as.character(Biostrings::translate(Biostrings::DNAStringSet(nt_seq), no.init.codon=T))])
  #Calculate amino acid hamming distance
  variant_dt[indel==F, Nham_aa := mapply(dimsum__hamming_distance, aa_seq, wt_AAseq)]
  #Reorder
  variant_dt <- variant_dt[,.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "Nham_nt", "Nham_aa", 
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel")]

  ### Indel variants
  ###########################

  if(!dimsum_meta[["indels"]]){
    variant_dt[indel==T, indel_discarded := T]
  }else{
    if(is.na(dimsum_meta[["indelLengths"]][1])){
      variant_dt[indel==T, indel_discarded := F]
    }else{
      variant_dt[indel==T, indel_discarded := !nchar(nt_seq) %in% dimsum_meta[["indelLengths"]]]
    }
  }

  #Dicarded indel variant count statistics
  nuc_indel_dict <- as.list(apply(
    variant_dt[indel_discarded==T,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]], 
    2, sum))

  #Retained indel variant count statistics
  nuc_indrt_dict <- as.list(apply(
    variant_dt[indel_discarded==F,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]], 
    2, sum))

  #Indel variants
  indel_variant_dt <- data.table::copy(variant_dt[indel_discarded==T,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel_discarded")])

  ### Update remaining variants
  ###########################

  #Remaining variants without discarded indels
  variant_dt <- variant_dt[indel==F | indel_discarded==F,.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "Nham_nt", "Nham_aa",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel")]
  #Indicate WT sequence
  variant_dt[nt_seq == wt_ntseq,WT := TRUE]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after indel variant filtering"), call. = FALSE)
  }

  #Check if WT sequence exists
  if(variant_dt[WT==T,.N]==0){
    variant_dt[, all_reads := rowSums(.SD > 0) == length(grep("_count$", names(variant_dt))),,.SDcols = grep("_count$", names(variant_dt))]
    variant_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("_s[0].*_count$"), names(variant_dt))]
    dimsum__status_message(paste0("WT variant not found. Did you mean to specify one of the following?\n"))
    print(variant_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
    stop(paste0("Cannot proceed with variant processing: WT variant not found"), call. = FALSE)
  }

  #Indicate STOPs
  variant_dt <- dimsum__identify_STOP_mutations(variant_dt)

  ### Rejected variants (mutated constant region)
  ###########################

  if(sum(!strsplit(dimsum_meta[["wildtypeSequenceCoded"]], "")[[1]] %in% c("A", "C", "G", "T"))!=0){
    variant_dt <- dimsum__remove_internal_constant_region(
      input_dt = variant_dt,
      wt_ccntseq = dimsum_meta[["wildtypeSequenceCoded"]])
  }else{
    variant_dt[, constant_region := T]
  }

  #Mutated constant region variant count statistics
  nuc_const_dict <- as.list(apply(
    variant_dt[constant_region==F,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]], 
    2, sum))

  #Constant region variants
  rejected_variant_dt <- data.table::copy(variant_dt[constant_region==F,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))], 
    "barcode_valid", "indel", "constant_region")])

  ### Update remaining variants
  ###########################

  #Remaining variants without mutated constant region (or with indels if retained)
  variant_dt <- variant_dt[constant_region==T | (indel==T & dimsum_meta[["indels"]]),]
  #Indicate STOPs
  variant_dt <- dimsum__identify_STOP_mutations(variant_dt)
  #Recalculate nucleotide hamming distance
  variant_dt[indel==F, Nham_nt := mapply(dimsum__hamming_distance, nt_seq, variant_dt[WT==T,nt_seq])]
  #Recalculate amino acid hamming distance
  variant_dt[indel==F, Nham_aa := mapply(dimsum__hamming_distance, aa_seq, variant_dt[WT==T,aa_seq])]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after mutated constant region filtering"), call. = FALSE)
  }

  ### Rejected variants (forbidden mutations)
  ###########################

  variant_dt <- dimsum__identify_permitted_mutations(
    dimsum_meta = dimsum_meta,
    input_dt = variant_dt)

  #Forbidden variant count statistics
  nuc_frbdn_dict <- as.list(apply(
    variant_dt[permitted==F,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]], 
    2, sum))

  #Forbidden variants
  rejected_variant_dt <- rbind(rejected_variant_dt, 
    data.table::copy(variant_dt[permitted==F,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted")]), fill = T)

  ### Update remaining variants
  ###########################

  #Remaining variants with permitted mutations (or with indels if retained)
  variant_dt <- variant_dt[permitted==T | (indel==T & dimsum_meta[["indels"]]),]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after forbidden mutation filtering"), call. = FALSE)
  }

  ### Rejected variants (too many substitutions)
  ###########################

  variant_dt[indel==F, too_many_substitutions := F]
  if(dimsum_meta[["sequenceType"]]=="coding"){
    variant_dt[indel==F & Nham_aa>dimsum_meta[["maxSubstitutions"]], too_many_substitutions := T]
  }else{
    variant_dt[indel==F & Nham_nt>dimsum_meta[["maxSubstitutions"]], too_many_substitutions := T]
  }

  #Too many substitution variant count statistics
  nuc_tmsub_dict <- as.list(apply(
    variant_dt[too_many_substitutions==T,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]], 
    2, sum))

  #Too many substitution variants
  rejected_variant_dt <- rbind(rejected_variant_dt, 
    data.table::copy(variant_dt[too_many_substitutions==T,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted", "too_many_substitutions")]), fill = T)

  ### Update remaining variants
  ###########################

  #Remaining variants without too many substitutions (or with indels if retained)
  variant_dt <- variant_dt[too_many_substitutions==F | (indel==T & dimsum_meta[["indels"]]),]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after too many substitutions filtering"), call. = FALSE)
  }

  ### Rejected variants (mixed substitutions)
  ###########################

  #Check if WT sequence exists
  if(variant_dt[WT==T,.N]==0){
    variant_dt[, all_reads := rowSums(.SD > 0) == length(grep("_count$", names(variant_dt))),,.SDcols = grep("_count$", names(variant_dt))]
    variant_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("_s[0].*_count$"), names(variant_dt))]
    dimsum__status_message(paste0("WT variant not found. Did you mean to specify one of the following?\n"))
    print(variant_dt[all_reads == T,][order(mean_count, decreasing = T)[1:5],.(nt_seq = toupper(nt_seq), all_reads, mean_count)])
    stop(paste0("Cannot proceed with variant processing: WT variant not found"), call. = FALSE)
  }

  #Add number of codons affected by mutations
  wt_ntseq <- variant_dt[WT==T,nt_seq]
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]
  variant_dt[indel==F, Nmut_codons := length(unique(ceiling(which(strsplit(nt_seq,"")[[1]] != wt_ntseq_split)/3))),nt_seq]

  #Identify mixed substitutions
  variant_dt[indel==F, mixed_substitutions := FALSE]
  if(dimsum_meta[["sequenceType"]]=="coding" & !dimsum_meta[["mixedSubstitutions"]]){
    variant_dt[indel==F & (Nmut_codons-Nham_aa) != 0 & Nham_aa != 0, mixed_substitutions := TRUE]
  }

  #Mixed substitution variant count statistics
  nuc_mxsub_dict <- as.list(apply(
    variant_dt[mixed_substitutions==T,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]], 
    2, sum))

  #Mixed substitution variants
  rejected_variant_dt <- rbind(rejected_variant_dt, 
    data.table::copy(variant_dt[mixed_substitutions==T,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted", "too_many_substitutions", "mixed_substitutions")]), fill = T)

  ### Remaining variants
  ###########################

  #Remaining variants without mixed substitutions (or with indels if retained)
  variant_dt <- variant_dt[mixed_substitutions==F | (indel==T & dimsum_meta[["indels"]]),.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "WT", "STOP", "STOP_readthrough", "Nham_nt", "Nham_aa", "Nmut_codons",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted", "too_many_substitutions", "mixed_substitutions")]

  #Check if variants remaining
  if(variant_dt[,.N]==0){
    stop(paste0("Cannot proceed with variant processing: No variants found after mixed substitutions filtering"), call. = FALSE)
  }

  #Filtered variant nucleotide distributions
  for(count_column in names(variant_dt)[grepl("_count$", names(variant_dt))]){
    nuc_subst_dict[[count_column]] <- tapply(
      unlist(variant_dt[,.SD,,.SDcols = count_column]), 
      variant_dt[, Nham_nt], 
      sum, na.rm = T)
  }

  #Filtered variant amino acid distributions
  for(count_column in names(variant_dt)[grepl("_count$", names(variant_dt))]){
    aa_subst_dict[[count_column]] <- tapply(
      unlist(variant_dt[,.SD,,.SDcols = count_column]), 
      variant_dt[, Nham_aa], 
      sum, na.rm = T)
  }

  ### Save variant statistics
  ###########################

  #Save mutation statistics dictionaries
  mutation_stats_dicts <- list(
    "aa_subst_dict" = aa_subst_dict,
    "nuc_subst_dict" = nuc_subst_dict, 
    "nuc_indrt_dict" = nuc_indrt_dict,
    "nuc_mxsub_dict" = nuc_mxsub_dict,
    "nuc_tmsub_dict" = nuc_tmsub_dict,
    "nuc_frbdn_dict" = nuc_frbdn_dict, 
    "nuc_const_dict" = nuc_const_dict, 
    "nuc_indel_dict" = nuc_indel_dict,
    "nuc_nbarc_dict" = nuc_nbarc_dict)
  save(mutation_stats_dicts, 
    file = file.path(dimsum_meta[["tmp_path"]], "mutation_stats_dicts.RData"),
    version = 2)

  dimsum__status_message("Done\n")

  return(list(
    "nobarcode_variants" = nobarcode_variant_dt,
    "indel_variants" = indel_variant_dt,
    "rejected_variants" = rejected_variant_dt,
    "retained_variants" = variant_dt))

}
