
#' dimsum__merge_fitness
#'
#' Calculate fitness and count-based error.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param doubles_dt doubles data.table (required)
#' @param singles_dt singles data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param fitness_outpath output path for saved objects (required)
#' @param report whether or not to generate fitness summary plots (default: TRUE)
#' @param report_outpath fitness report output path
#'
#' @return Nothing
#' @export
#' @import data.table
dimsum__merge_fitness <- function(
  dimsum_meta,
  input_dt,
  doubles_dt,
  singles_dt,
  all_reps,
  fitness_outpath,
  report = TRUE,
  report_outpath = NULL
  ){

  dimsum__status_message("Merging fitness estimates from biological replicates...\n")

  aa_obj <- Biostrings::AAString("GAVLMIFYWKRHDESTCNQP")
  aa_list <- Biostrings::AMINO_ACID_CODE[strsplit(as.character(aa_obj), NULL)[[1]]]
  aa_list["*"] <- "*"

  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")

  #### all variants
  fitness_rx <- input_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]"),colnames(input_dt))]
  sigma_rx <- input_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]"),colnames(input_dt))]
  input_dt[,fitness := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  input_dt[,sigma := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]

  #### singles
  if(nrow(singles_dt)!=0){
    fitness_rx <- singles_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]"),colnames(singles_dt))]
    sigma_rx <- singles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]"),colnames(singles_dt))]
    singles_dt[,fitness := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
    singles_dt[,sigma := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  }

  #### doubles
  #uncorrected fitness
  if(nrow(doubles_dt)!=0){
    fitness_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]_uncorr"),colnames(doubles_dt))]
    sigma_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_uncorr"),colnames(doubles_dt))]
    doubles_dt[,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
    doubles_dt[,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  }

  #conditioned fitness
  if(dimsum_meta[["bayesianDoubleFitness"]] & nrow(doubles_dt)!=0){
    fitness_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]_cond"),colnames(doubles_dt))]
    sigma_rx <- doubles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_cond"),colnames(doubles_dt))]
    doubles_dt[,fitness_cond := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
    doubles_dt[,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  }

  dimsum__status_message("Done\n")

  ### Add growth rates if cell_density and selection_time (in hours) supplied
  ###########################

  #Check that selection_time supplied for all output samples
  if(sum(is.na(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,"selection_time"]))==0){
    #Check that cell_density supplied for all samples
    if(sum(is.na(dimsum_meta[["exp_design"]][,"cell_density"]))==0){

      dimsum__status_message("Inferring growth rates...\n")

      input_dt <- dimsum__infer_growth_rates(
        dimsum_meta = dimsum_meta,
        input_dt = input_dt,
        all_reps = all_reps
        )

      dimsum__status_message("Done\n")

    }
  }

  ### Output replicate data files
  ###########################

  dimsum__status_message("Saving fitness estimates...\n")

  #Rename columns
  names(input_dt)[names(input_dt)=="merge_seq"] <- "nt_seq"
  names(singles_dt)[names(singles_dt)=="merge_seq"] <- "nt_seq"
  names(doubles_dt)[names(doubles_dt)=="merge_seq"] <- "nt_seq"

  #Reformat columns
  if(dimsum_meta[["sequenceType"]]=="coding"){
    #Set nucleotide sequence to NA when fitness and error aggregated at AA level
    if(dimsum_meta[["mixedSubstitutions"]]){
      input_dt[, nt_seq := NA]
      if(nrow(singles_dt)!=0){singles_dt[, nt_seq := NA]}
    }else{
      input_dt[Nham_aa>0 | indel==T, nt_seq := NA]
      if(nrow(singles_dt)!=0){singles_dt[Nham_aa>0, nt_seq := NA]}
    }
    if(nrow(doubles_dt)!=0){doubles_dt[, nt_seq := NA]}
    #Remove unnecessary columns
    unnecessary_cols <- c(
      "var_fitness",
      "avg_sigma",
      "isNA",
      "avg_sigma_fit",
      "merge_seq",
      "counts_for_bins",
      "bin_count")
    input_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% unnecessary_cols]]
    singles_dt <- singles_dt[,.SD,,.SDcols = names(singles_dt)[!names(singles_dt) %in% unnecessary_cols]]
    doubles_dt <- doubles_dt[,.SD,,.SDcols = names(doubles_dt)[!names(doubles_dt) %in% unnecessary_cols]]
  }else{
    #Remove unnecessary columns
    unnecessary_cols <- c(
      "aa_seq",
      "Nham_aa",
      "Nmut_codons",
      "STOP",
      "STOP_readthrough",
      "var_fitness",
      "avg_sigma",
      "isNA",
      "avg_sigma_fit",
      "merge_seq",
      "counts_for_bins",
      "bin_count")
    input_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[!names(input_dt) %in% unnecessary_cols]]
    singles_dt <- singles_dt[,.SD,,.SDcols = names(singles_dt)[!names(singles_dt) %in% unnecessary_cols]]
    doubles_dt <- doubles_dt[,.SD,,.SDcols = names(doubles_dt)[!names(doubles_dt) %in% unnecessary_cols]]
    #Rename WT_AA columns to WT_nt
    names(input_dt)[grep("^WT_AA", names(input_dt))] <- gsub("WT_AA", "WT_nt", names(input_dt)[grep("^WT_AA", names(input_dt))])
    names(singles_dt)[grep("^WT_AA", names(singles_dt))] <- gsub("WT_AA", "WT_nt", names(singles_dt)[grep("^WT_AA", names(singles_dt))])
    names(doubles_dt)[grep("^WT_AA", names(doubles_dt))] <- gsub("WT_AA", "WT_nt", names(doubles_dt)[grep("^WT_AA", names(doubles_dt))])
  }

  #Rename objects
  all_variants <- input_dt
  wildtype <- all_variants[WT==T,,,.SDcols = names(all_variants)[!grepl("^count_", names(all_variants))]]
  doubles <- doubles_dt

  ### Output plain text files
  ###########################

  ##### finalize data.tables
  silent <- data.table()
  singles <- data.table()
  singles_mavedb <- data.table()
  if(nrow(singles_dt)!=0){
    if(dimsum_meta[["sequenceType"]]=="coding"){
      silent <- singles_dt[Nham_aa==0,.(Pos,WT_AA,Mut,nt_seq,aa_seq,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,mean_count,fitness,sigma)]
      singles <- singles_dt[Nham_aa==1,.(Pos,WT_AA,Mut,nt_seq,aa_seq,Nham_nt,Nham_aa,Nmut_codons,STOP,STOP_readthrough,mean_count,fitness,sigma)]
      singles_mavedb <- singles[,.(hgvs_pro = NA, score = fitness, se = sigma)]
      singles_mavedb[, hgvs_pro := paste0("p.", aa_list[singles[,WT_AA]], singles[,Pos], aa_list[singles[,Mut]])]
    }else{
      singles <- singles_dt[,.(Pos,WT_nt,Mut,nt_seq,Nham_nt,mean_count,fitness,sigma)]
      singles_mavedb <- singles[,.(hgvs_nt = NA, score = fitness, se = sigma)]
      singles_mavedb[, hgvs_nt := paste0("n.", singles[,Pos], toupper(singles[,WT_nt]), ">", toupper(singles[,Mut]))]
    }
  }

  #for doubles #add single mutant fitness/sigma values to double mutant table
  if(nrow(doubles_dt)!=0 & nrow(singles)!=0){
    doubles_dt[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
    doubles_dt[,sigma1 := singles[Pos == Pos1 & Mut == Mut1,sigma],.(Pos1,Mut1)]
    doubles_dt[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
    doubles_dt[,sigma2 := singles[Pos == Pos2 & Mut == Mut2,sigma],.(Pos2,Mut2)]
  }

  if(dimsum_meta[["sequenceType"]]=="coding"){
    retained_cols <- c("Pos1","Pos2","WT_AA1","WT_AA2","Mut1","Mut2","nt_seq","aa_seq","Nham_nt","Nham_aa","Nmut_codons","STOP","STOP_readthrough","mean_count",
      "fitness1","sigma1","fitness2","sigma2",
      "fitness_uncorr","sigma_uncorr",
      "fitness_cond","sigma_cond")
    doubles <- doubles_dt[,.SD,,.SDcols = retained_cols[retained_cols %in% names(doubles_dt)]]

    #write data to files
    write.table(x = wildtype, file = file.path(fitness_outpath, "fitness_wildtype.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = silent, file = file.path(fitness_outpath, "fitness_silent.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = singles, file = file.path(fitness_outpath, "fitness_singles.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = doubles, file = file.path(fitness_outpath, "fitness_doubles.txt"),
                quote = F,row.names = F, col.names = T)

  }else{
    retained_cols <- c("Pos1","Pos2","WT_nt1","WT_nt2","Mut1","Mut2","nt_seq","Nham_nt","mean_count",
      "fitness1","sigma1","fitness2","sigma2",
      "fitness_uncorr","sigma_uncorr",
      "fitness_cond","sigma_cond")
    doubles <- doubles_dt[,.SD,,.SDcols = retained_cols[retained_cols %in% names(doubles_dt)]]

    #write data to files
    write.table(x = wildtype, file = file.path(fitness_outpath, "fitness_wildtype.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = singles, file = file.path(fitness_outpath, "fitness_singles.txt"),
                quote = F,row.names = F, col.names = T)
    write.table(x = doubles, file = file.path(fitness_outpath, "fitness_doubles.txt"),
                quote = F,row.names = F, col.names = T)
  }

  #Write MaveDB formatted singles
  write.table(x = singles_mavedb, file = file.path(fitness_outpath, "fitness_singles_MaveDB.csv"),
              quote = F,row.names = F, col.names = T, sep = ",")

  ### Output RData file
  ###########################

  if(dimsum_meta[["sequenceType"]]=="coding"){
    save(all_variants, wildtype, silent, singles, doubles, 
      file = file.path(fitness_outpath, paste0(dimsum_meta[["projectName"]], '_fitness_replicates.RData')),
      version = 2)
  }else{
    save(all_variants, wildtype, singles, doubles, 
      file = file.path(fitness_outpath, paste0(dimsum_meta[["projectName"]], '_fitness_replicates.RData')),
      version = 2)
  }

  dimsum__status_message("Done\n")

}
