
#' dimsum_stage_counts_to_fitness
#'
#' Estimate fitness of singles, doubles and silent AA substitution variants from DiMSum formatted output counts.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param fitness_outpath output path for saved objects (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum_stage_counts_to_fitness <- function(
  dimsum_meta,
  fitness_outpath,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){

  #Whether or not to execute the system command
  this_stage <- 5
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #This stage after stopStage
  if(dimsum_meta[["stopStage"]] < this_stage){
    return(dimsum_meta)
  }
  
  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Fitness path
  dimsum_meta[["fitness_path"]] <- file.path(fitness_outpath)

  #Do nothing if analysis not executed
  if(!execute){
    return(dimsum_meta)
  }

  #Create/overwrite fitness directory (if executed)
  fitness_outpath <- gsub("/$", "", fitness_outpath)
  dimsum__create_dir(fitness_outpath, execute = execute, message = "DiMSum STAGE 5 (STEAM): ANALYSE VARIANT COUNTS", overwrite_dir = FALSE) 

  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = dimsum_meta[["variant_data_merge_path"]],
    stage_number = this_stage,
    execute = execute)

  ### Load variant data
  ###########################
    
  #load variant data from RData file
  e1 <- new.env() 
  load(dimsum_meta[["variant_data_merge_path"]])
  all_data <- get('variant_data_merge', e1)
  dimsum__check_variants(dimsum_meta = dimsum_meta, input_dt = all_data)
  rm(e1)

  #WT nucleotide sequences
  wt_ntseq <- all_data[WT==T,nt_seq]

  #WT AA sequences
  wt_AAseq <- all_data[WT==T,aa_seq]

  #Determine replicates to retain
  all_reps <- unique(dimsum_meta[["exp_design"]][,"experiment"])
  if(dimsum_meta[["retainedReplicates"]]!="all"){
    all_reps <- unique(as.integer(strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]]))
  }

  #Check if only one experimental replicate exists (requirement for normalisation and error model fit)
  if(length(unique(dimsum_meta[["exp_design"]][,"experiment_replicate"]))<2){
    dimsum_meta[["fitnessNormalise"]] <- FALSE
    dimsum_meta[["fitnessErrorModel"]] <- FALSE
  }

  ### Filter out low count nucleotide variants
  ###########################

  nf_data <- dimsum__filter_nuc_variants(
    dimsum_meta = dimsum_meta,
    input_dt = all_data,
    all_reps = all_reps)

  ### Identify variants to use for error modelling
  ###########################

  nf_data[, error_model := T]

  ### Aggregate counts at the AA level (if coding sequence and "mixedSubstitutions"==T)
  ###########################

  #For coding sequences aggregate variant counts at the AA level
  if(dimsum_meta[["sequenceType"]]=="coding" & dimsum_meta[["mixedSubstitutions"]]){
    nf_data <- dimsum__aggregate_AA_variants(
      dimsum_meta = dimsum_meta,      
      input_dt = nf_data,
      all_reps = all_reps)
  }else{
    nf_data[,merge_seq := nt_seq]
  }

  ### Aggregate counts for biological output replicates
  ###########################

  dimsum__status_message("Aggregating counts for biological output replicates...\n")

  #Add up counts for biological output reps
  for (E in all_reps) {
    idx <- names(nf_data)[grep(names(nf_data),pattern = paste0("e",E,"_s1_b"))]
    nf_data[!rowSums(is.na(nf_data[,.SD,.SDcols = idx]))==length(idx),paste0("count_e",E,"_s1") := rowSums(.SD, na.rm = T),,.SDcols = idx]
    names(nf_data)[grep(names(nf_data),pattern = paste0("e",E,"_s0_b"))] <- paste0("count_e",E,"_s0")
  }
  nf_data <- nf_data[,.SD,,.SDcols = c(
    "merge_seq","nt_seq","aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","indel","STOP","STOP_readthrough","error_model",names(nf_data)[grep(names(nf_data),pattern="^count")])]

  dimsum__status_message("Done\n")

  ### Fit error model
  ###########################

  #Fit error model
  model_result <- dimsum__error_model(
    dimsum_meta = dimsum_meta,
    input_dt = data.table::copy(nf_data[error_model==T]),
    all_reps = all_reps,
    report_outpath = report_outpath)

  ### Calculate fitness and count-based error (and remove variants without fitness estimates in any replicate)
  ###########################

  nff_data <- dimsum__calculate_fitness(
    dimsum_meta = dimsum_meta,
    input_dt = nf_data,
    all_reps = all_reps,
    error_model_dt = model_result[["error_model"]],
    norm_model_dt = model_result[["norm_model"]])

  ### Merge fitness and error at the AA level (if coding sequence and "mixedSubstitutions"==F)
  ###########################

  #Aggregate nonsynonymous variant fitness and error at the AA level
  if(dimsum_meta[["sequenceType"]]=="coding" & !dimsum_meta[["mixedSubstitutions"]]){
    nff_data <- dimsum__aggregate_AA_variants_fitness(
      dimsum_meta = dimsum_meta,      
      input_dt = nff_data,
      all_reps = all_reps)
  }else{
    if(dimsum_meta[["sequenceType"]]!="coding"){nff_data[,merge_seq := nt_seq,nt_seq]}
    nff_data <- nff_data[,.SD,merge_seq,.SDcols = c(
      "aa_seq","Nham_nt","Nham_aa",
      "Nmut_codons","WT","indel","STOP","STOP_readthrough", "error_model",
      names(nff_data)[grep(names(nff_data),pattern="^count")],
      "mean_count",
      names(nff_data)[grep(names(nff_data),pattern="^fitness|sigma")])]
  }

  ### Wild type
  ###########################

  wildtype <- nff_data[WT==TRUE & error_model==TRUE,]

  ### Identify position and identity of single AA/NT mutations (and all silent mutants in the case of AA mutations)
  ###########################

  if(dimsum_meta[["sequenceType"]]=="coding"){
    singles <- dimsum__identify_single_aa_mutations(
      input_dt = nff_data[error_model==T],
      wt_AAseq = wt_AAseq)
  }else{
    singles <- dimsum__identify_single_nt_mutations(
      input_dt = nff_data[error_model==T],
      wt_ntseq = wt_ntseq)    
  }

  ### Identify position and identity of double AA/NT mutations
  ###########################

  if(dimsum_meta[["sequenceType"]]=="coding"){
    doubles <- dimsum__identify_double_aa_mutations(
      input_dt = nff_data[error_model==T],
      singles_dt = singles,
      wt_AAseq = wt_AAseq)
  }else{
    doubles <- dimsum__identify_double_nt_mutations(
      input_dt = nff_data[error_model==T],
      singles_dt = singles,
      wt_ntseq = wt_ntseq)
  }

  ### Bayesian framework for fitness estimation for double mutants
  ###########################

  #Bin mean counts for replicate 1
  # doubles[,counts_for_bins := .SD[[1]],,.SDcols = paste0("count_e",all_reps[1],"_s0")]
  # doubles[,bin_count := findInterval(log10(counts_for_bins),seq(0.5,4,0.25))]
  # doubles[,.(.N,mean(counts_for_bins)),bin_count][order(bin_count)]

  if(!dimsum_meta[["bayesianDoubleFitness"]]){
    dimsum__status_message("Skipping Bayesian double mutant fitness estimation\n")
  }else{
    doubles <- dimsum__bayesian_double_fitness(
      dimsum_meta = dimsum_meta,
      doubles_dt = doubles,
      singles_dt = singles,
      wt_dt = wildtype,
      all_reps = all_reps,
      report_outpath = report_outpath)
  }

  ### Normalise fitness and sigma for differences in number of generations (between biological replicates)
  ###########################

  #Check that number of generations supplied for all output samples
  if(sum(is.na(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,"generations"]))!=0){

    dimsum__status_message("Skipping normalisation by number of generations\n")

  }else{

    dimsum__status_message("Normalising fitness and error for differences in number of generations...\n")

    nff_data <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = nff_data,
      all_reps = all_reps,
      fitness_suffix="_uncorr"
      )
    singles <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = singles,
      all_reps = all_reps
      )
    doubles <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = doubles,
      all_reps = all_reps,
      fitness_suffix="_uncorr"
      )
    if(dimsum_meta[["bayesianDoubleFitness"]]){
      doubles <- dimsum__normalise_fitness(
        dimsum_meta = dimsum_meta,
        input_dt = doubles,
        all_reps = all_reps,
        fitness_suffix="_cond"
        )
    }

    dimsum__status_message("Done\n")

  }

  ### Merge fitness values and save fitness data
  ###########################

  dimsum__merge_fitness(
    dimsum_meta = dimsum_meta,
    input_dt = nff_data,
    doubles_dt = doubles,
    singles_dt = singles,
    all_reps = all_reps,
    fitness_outpath = fitness_outpath,
    report = TRUE,
    report_outpath = report_outpath)

  #Delete files when last stage complete
  if(!dimsum_meta[["retainIntermediateFiles"]]){
    if(dimsum_meta[["stopStage"]]==this_stage){
      if(!is.null(dimsum_meta[["deleteIntermediateFiles"]])){suppressWarnings(temp_out <- file.remove(dimsum_meta[["deleteIntermediateFiles"]]))}
      if(!is.null(dimsum_meta[["touchIntermediateFiles"]])){suppressWarnings(temp_out <- file.create(dimsum_meta[["touchIntermediateFiles"]]))}
    }
  }

  return(dimsum_meta)
}
