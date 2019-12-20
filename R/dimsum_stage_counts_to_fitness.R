
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
  this_stage <- 7
  execute <- (dimsum_meta[["startStage"]] <= this_stage & (dimsum_meta[["stopStage"]] == 0 | dimsum_meta[["stopStage"]] >= this_stage))
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
  dimsum__create_dir(fitness_outpath, execute = execute, message = "DiMSum STAGE 7: CALCULATE FITNESS", overwrite_dir = FALSE) 

  ### Load variant data
  ###########################

  #load variant data from RData file
  e1 <- new.env() 
  load(dimsum_meta[["variant_data_merge_path"]])
  all_data <- get('variant_data_merge', e1)
  rm(e1)

  #WT nucleotide sequences
  wt_ntseq <- all_data[WT==T,nt_seq]

  #WT AA sequences
  wt_AAseq <- all_data[WT==T,aa_seq]

  #Sample names
  input_samples <- names(all_data)[grep("_e.*_s0_b.*_count$", names(all_data))]
  output_samples <- names(all_data)[grep("_e.*_s1_b.*_count$", names(all_data))]

  #Determine replicates to retain
  all_reps <- unique(dimsum_meta[["exp_design"]][,"experiment"])
  if(dimsum_meta[["retainedReplicates"]]!="all"){
    all_reps <- unique(as.integer(strsplit(dimsum_meta[["retainedReplicates"]], ",")[[1]]))
  }

  #Bayesian double mutant fitness estimates
  bayesian_double_fitness <- dimsum_meta[["bayesianDoubleFitness"]]

  ### Filter out low count nucleotide variants
  ###########################

  nf_data <- dimsum__filter_nuc_variants(
    dimsum_meta = dimsum_meta,
    input_dt = all_data,
    all_reps = all_reps)

  ### Fit error model
  ###########################

  #Fit error model (using variants with less than specified number of mutations)
  model_result <- dimsum__error_model(
    dimsum_meta = dimsum_meta,
    input_dt = data.table::copy(nf_data),
    all_reps = all_reps,
    report_outpath = report_outpath)

  ### Aggregate counts from variants that are identical at the AA level and without synonymous mutations (if coding sequence)
  ###########################

  #Add number of codons affected by mutations
  wt_ntseq_split <- strsplit(wt_ntseq,"")[[1]]
  nf_data[,Nmut_codons := length(unique(ceiling(which(strsplit(nt_seq,"")[[1]] != wt_ntseq_split)/3))),nt_seq]
  #For coding sequences retain only either purely nonsynonymous mutations or purely silent mutations
  if(dimsum_meta[["sequenceType"]]=="coding" & !dimsum_meta[["mixedSubstitutions"]]){
    nf_data <- nf_data[(Nmut_codons-Nham_aa) == 0 | Nham_aa == 0,]
  }

  if(dimsum_meta[["sequenceType"]]=="coding"){
    nf_data_syn <- dimsum__aggregate_AA_variants(
      input_dt = nf_data)
  }else{
    nf_data[,merge_seq := nt_seq,nt_seq]
    nf_data_syn <- nf_data[,.SD,merge_seq,.SDcols = c(
      "aa_seq","Nham_nt","Nham_aa",
      "Nmut_codons","WT","STOP","STOP_readthrough",names(nf_data)[grep(names(nf_data),pattern="_count$")])]
  }

  ### Aggregate counts for biological output replicates
  ###########################

  message("Aggregating counts for biological output replicates...")

  #Add up counts for biological output reps
  for (E in all_reps) {
    idx <- names(nf_data_syn)[grep(names(nf_data_syn),pattern = paste0("e",E,"_s1_b"))]
    nf_data_syn[,paste0("count_e",E,"_s1") := rowSums(.SD),,.SDcols = idx]
    names(nf_data_syn)[grep(names(nf_data_syn),pattern = paste0("e",E,"_s0_b"))] <- paste0("count_e",E,"_s0")
  }
  nf_data_syn <- unique(nf_data_syn[,.SD,merge_seq,.SDcols = c(
    "aa_seq","Nham_nt","Nham_aa",
    "Nmut_codons","WT","STOP","STOP_readthrough",names(nf_data_syn)[grep(names(nf_data_syn),pattern="^count")])])

  message("Done")

  ### Calculate fitness and count-based error (and remove variants without fitness estimates in any replicate)
  ###########################

  nf_data_syn <- dimsum__calculate_fitness(
    dimsum_meta = dimsum_meta,
    input_dt = nf_data_syn,
    all_reps = all_reps,
    error_model_dt = model_result[["error_model"]],
    norm_model_dt = model_result[["norm_model"]])

  ### Wild type
  ###########################

  wildtype <- nf_data_syn[WT==TRUE,]

  ### Identify position and identity of single AA/NT mutations (and all silent mutants in the case of AA mutations)
  ###########################

  if(dimsum_meta[["sequenceType"]]=="coding"){
    singles_silent <- dimsum__identify_single_aa_mutations(
      input_dt = nf_data_syn,
      wt_AAseq = wt_AAseq)
  }else{
    singles_silent <- dimsum__identify_single_nt_mutations(
      input_dt = nf_data_syn,
      wt_ntseq = wt_ntseq)    
  }

  ### Identify position and identity of double AA/NT mutations
  ###########################

  if(dimsum_meta[["sequenceType"]]=="coding"){
    doubles <- dimsum__identify_double_aa_mutations(
      input_dt = nf_data_syn,
      singles_dt = singles_silent,
      wt_AAseq = wt_AAseq)
  }else{
    doubles <- dimsum__identify_double_nt_mutations(
      input_dt = nf_data_syn,
      singles_dt = singles_silent,
      wt_ntseq = wt_ntseq)
  }

  ### Bayesian framework for fitness estimation for double mutants
  ###########################

  #Bin mean counts for replicate 1
  doubles[,counts_for_bins := .SD[[1]],,.SDcols = paste0("count_e",all_reps[1],"_s0")]
  doubles[,bin_count := findInterval(log10(counts_for_bins),seq(0.5,4,0.25))]
  # doubles[,.(.N,mean(counts_for_bins)),bin_count][order(bin_count)]

  if(!bayesian_double_fitness){
    message("Skipping Bayesian double mutant fitness estimation")
  }else{
    doubles <- dimsum__bayesian_double_fitness(
      dimsum_meta = dimsum_meta,
      doubles_dt = doubles,
      singles_dt = singles_silent,
      wt_dt = wildtype,
      all_reps = all_reps,
      report_outpath = report_outpath)
  }

  ### Normalise fitness and sigma for differences in number of generations (between biological replicates)
  ###########################

  #Check that number of generations supplied for all output samples
  if(sum(is.na(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,"generations"]))!=0){

    message("Skipping normalisation by number of generations")

  }else{

    message("Normalising fitness and error for differences in number of generations...")

    nf_data_syn <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = nf_data_syn,
      all_reps = all_reps,
      fitness_suffix="_uncorr"
      )
    singles_silent <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = singles_silent,
      all_reps = all_reps
      )
    doubles <- dimsum__normalise_fitness(
      dimsum_meta = dimsum_meta,
      input_dt = doubles,
      all_reps = all_reps,
      fitness_suffix="_uncorr"
      )
    if(bayesian_double_fitness){
      doubles <- dimsum__normalise_fitness(
        dimsum_meta = dimsum_meta,
        input_dt = doubles,
        all_reps = all_reps,
        fitness_suffix="_cond"
        )
    }

    message("Done")

  }

  ### Merge fitness values and save fitness data
  ###########################

  dimsum__merge_fitness(
    dimsum_meta = dimsum_meta,
    input_dt = nf_data_syn,
    doubles_dt = doubles,
    singles_dt = singles_silent,
    wt_dt = wildtype,
    all_reps = all_reps,
    fitness_outpath = fitness_outpath,
    report = TRUE,
    report_outpath = report_outpath)

  #Delete files when last stage complete
  if(!dimsum_meta[["retainIntermediateFiles"]]){
    if(dimsum_meta[["stopStage"]]==this_stage | dimsum_meta[["stopStage"]]==0){
      temp_out <- mapply(system, dimsum_meta[["deleteIntermediateFiles"]], MoreArgs = list(ignore.stdout = T, ignore.stderr = T))
    }
  }

  return(dimsum_meta)
}
