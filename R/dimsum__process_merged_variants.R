
#' dimsum__process_merged_variants
#'
#' Process nucleotide variants to remove indels, internal constant region variants, not permitted sequences and variants with too many substitions.
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

  message("Processing merged variants...")

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

  #Define indel variants as those with length different from WT
  variant_dt[, indel := F]
  variant_dt[nchar(nt_seq)!=nchar(wt_ntseq), indel := T]
  #Calculate nucleotide hamming distance
  variant_dt[indel==F, Nham_nt := mapply(dimsum__hamming_distance, nt_seq, wt_ntseq)]

  #Determine amino acid sequence
  suppressWarnings(variant_dt[, aa_seq := as.character(Biostrings::translate(Biostrings::DNAStringSet(nt_seq)))])
  #Calculate amino acid hamming distance
  variant_dt[indel==F, Nham_aa := mapply(dimsum__hamming_distance, aa_seq, wt_AAseq)]
  #Reorder
  variant_dt <- variant_dt[,.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "Nham_nt", "Nham_aa", 
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel")]

  ### Indel variants
  ###########################

  #Indel variant count statistics
  nuc_indel_dict <- as.list(apply(
    variant_dt[indel==T,.SD,,.SDcols = names(variant_dt)[grepl("_count$", names(variant_dt))]], 
    2, sum))

  #Indel variants
  indel_variant_dt <- data.table::copy(variant_dt[indel==T,.SD,,.SDcols = c(
    "nt_seq", "aa_seq",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel")])

  ### Update remaining variants
  ###########################

  #Remaining variants without indels
  variant_dt <- variant_dt[indel==F,.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "Nham_nt", "Nham_aa",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel")]
  #Indicate WT sequence
  variant_dt[nt_seq == wt_ntseq,WT := TRUE]

  #Check if WT sequence exists
  if(variant_dt[WT==T,.N]==0){
    variant_dt[, all_reads := rowSums(.SD > 0) == length(grep("_count$", names(variant_dt))),,.SDcols = grep("_count$", names(variant_dt))]
    variant_dt[, mean_count := rowMeans(.SD),,.SDcols = grep(paste0("_s[0].*_count$"), names(variant_dt))]
    message(paste0("WT variant not found. Did you mean to specify one of the following?"))
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

  #Remaining variants without mutated constant region
  variant_dt <- variant_dt[constant_region==T,]
  #Indicate STOPs
  variant_dt <- dimsum__identify_STOP_mutations(variant_dt)
  #Recalculate nucleotide hamming distance
  variant_dt[, Nham_nt := mapply(dimsum__hamming_distance, nt_seq, variant_dt[WT==T,nt_seq])]
  #Recalculate amino acid hamming distance
  variant_dt[, Nham_aa := mapply(dimsum__hamming_distance, aa_seq, variant_dt[WT==T,aa_seq])]

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

  #Remaining variants with permitted mutations
  variant_dt <- variant_dt[permitted==T,]

  ### Rejected variants (too many substitutions)
  ###########################

  variant_dt[, too_many_substitutions := F]
  if(dimsum_meta[["sequenceType"]]=="coding"){
    variant_dt[Nham_aa>dimsum_meta[["maxSubstitutions"]], too_many_substitutions := T]
  }else{
    variant_dt[Nham_nt>dimsum_meta[["maxSubstitutions"]], too_many_substitutions := T]
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

  ### Remaining variants
  ###########################

  #Remaining variants without too many substition variants
  variant_dt <- variant_dt[too_many_substitutions==F,.SD,,.SDcols = c(
    "nt_seq", "aa_seq", "WT", "STOP", "STOP_readthrough", "Nham_nt", "Nham_aa",
    names(variant_dt)[grepl("_count$", names(variant_dt))],
    "barcode_valid", "indel", "constant_region", "permitted", "too_many_substitutions")]

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
    "nuc_tmsub_dict" = nuc_tmsub_dict,
    "nuc_frbdn_dict" = nuc_frbdn_dict, 
    "nuc_const_dict" = nuc_const_dict, 
    "nuc_indel_dict" = nuc_indel_dict,
    "nuc_nbarc_dict" = nuc_nbarc_dict)
  save(mutation_stats_dicts, 
    file = file.path(dimsum_meta[["tmp_path"]], "mutation_stats_dicts.RData"),
    version = 2)

  message("Done")

  return(list(
    "nobarcode_variants" = nobarcode_variant_dt,
    "indel_variants" = indel_variant_dt,
    "rejected_variants" = rejected_variant_dt,
    "retained_variants" = variant_dt))

}
