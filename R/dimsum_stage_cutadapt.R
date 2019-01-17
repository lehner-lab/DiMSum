
#' dimsum_stage_cutadapt
#'
#' Run cutadapt on all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param cutadapt_outpath cutadapt output path (required)
#' @param execute whether or not to execute the system command (default: TRUE)
#' @param report whether or not to generate cutadapt summary plots (default: TRUE)
#' @param report_outpath cutadapt report output path
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_cutadapt <- function(
  dimsum_meta,
  cutadapt_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){
  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}
  #Create/overwrite cutadapt directory (if executed)
  cutadapt_outpath <- gsub("/$", "", cutadapt_outpath)
  dimsum__create_dir(cutadapt_outpath, execute = execute, message = "DiMSum STAGE 3: TRIM CONSTANT REGIONS")  
  #Trim FASTQ file pairs
  message("Trimming FASTQ files with cutadapt:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]][,"pair_directory"], c(dimsum_meta[['exp_design']][,"pair1"], dimsum_meta[['exp_design']][,"pair2"]))
  print(all_fastq)
  message("Processing...")
  for(i in 1:dim(dimsum_meta[['exp_design']])[1]){
    #TODO: cutadapt binary path specifiable on commandline?
    message(paste0("\t", dimsum_meta[['exp_design']][i,c('pair1', 'pair2')]))
    #If not trans library: convert to linked adapters if 3' constant region expected to be sequenced
    if(!dimsum_meta[["transLibrary"]]){
      num_cut5f <- ifelse(is.na(dimsum_meta[['exp_design']][i,"cutadaptCut5First"]), 0, dimsum_meta[['exp_design']][i,"cutadaptCut5First"])
      num_cut5s <- ifelse(is.na(dimsum_meta[['exp_design']][i,"cutadaptCut5Second"]), 0, dimsum_meta[['exp_design']][i,"cutadaptCut5Second"])
      num_cut3f <- ifelse(is.na(dimsum_meta[['exp_design']][i,"cutadaptCut3First"]), 0, dimsum_meta[['exp_design']][i,"cutadaptCut3First"])
      num_cut3s <- ifelse(is.na(dimsum_meta[['exp_design']][i,"cutadaptCut3Second"]), 0, dimsum_meta[['exp_design']][i,"cutadaptCut3Second"])
      #Read 1
      if( (dimsum_meta[['exp_design']][i,"pair1_length"]-num_cut5f-num_cut3f) > (nchar(dimsum_meta[['exp_design']][i,"cutadapt5First"]) + nchar(dimsum_meta[['wildtypeSequence']])) ){
        dimsum_meta[['exp_design']][i,"cutadapt5First"] <- paste0(dimsum_meta[['exp_design']][i,"cutadapt5First"], "...", dimsum_meta[['exp_design']][i,"cutadapt3First"])
        dimsum_meta[['exp_design']][i,"cutadapt3First"] <- NA
      }
      #Read2 (use read1 length; read1 lengths can be variable due to inconsistent barcode trimming with cutadapt)
      if( (dimsum_meta[['exp_design']][i,"pair1_length"]-num_cut5s-num_cut3s) > (nchar(dimsum_meta[['exp_design']][i,"cutadapt5Second"]) + nchar(dimsum_meta[['wildtypeSequence']])) ){
        dimsum_meta[['exp_design']][i,"cutadapt5Second"] <- paste0(dimsum_meta[['exp_design']][i,"cutadapt5Second"], "...", dimsum_meta[['exp_design']][i,"cutadapt3Second"])
        dimsum_meta[['exp_design']][i,"cutadapt3Second"] <- NA
      }
    }
    #Check if this system command should be executed
    if(execute){
      #Options for removing constant regions from beginning or end of either read in pair
      temp_options = paste0(' -g ', dimsum_meta[['exp_design']][i,"cutadapt5First"], ' -G ', dimsum_meta[['exp_design']][i,"cutadapt5Second"])
      if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt3First"]) ){temp_options = paste0(temp_options, " -a ", dimsum_meta[['exp_design']][i,"cutadapt3First"])}
      if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt3Second"]) ){temp_options = paste0(temp_options, " -A ", dimsum_meta[['exp_design']][i,"cutadapt3Second"])}
      #Options for swapping read1 and read2
      temp_options_swap = paste0(' -g forward=', dimsum_meta[['exp_design']][i,"cutadapt5First"], ' -g reverse=', dimsum_meta[['exp_design']][i,"cutadapt5Second"])
      #Options for removing a fixed number of bases from beginning or end of either read in pair
      temp_cut_options = ''
      if( !is.na(dimsum_meta[['exp_design']][i,"cutadaptCut5First"]) ){temp_cut_options = paste0(temp_cut_options, " -u ", dimsum_meta[['exp_design']][i,"cutadaptCut5First"])}
      if( !is.na(dimsum_meta[['exp_design']][i,"cutadaptCut3First"]) ){temp_cut_options = paste0(temp_cut_options, " -u ", -dimsum_meta[['exp_design']][i,"cutadaptCut3First"])}
      if( !is.na(dimsum_meta[['exp_design']][i,"cutadaptCut5Second"]) ){temp_cut_options = paste0(temp_cut_options, " -U ", dimsum_meta[['exp_design']][i,"cutadaptCut5Second"])}
      if( !is.na(dimsum_meta[['exp_design']][i,"cutadaptCut3Second"]) ){temp_cut_options = paste0(temp_cut_options, " -U ", -dimsum_meta[['exp_design']][i,"cutadaptCut3Second"])}
      #Not stranded library
      if( !dimsum_meta[["stranded"]] ){
        #Swap reads in fastq files according to adapter presence
        temp_out = system(paste0(
          "cutadapt",
          temp_options_swap,
          temp_cut_options,
          " --no-trim ",
          " --minimum-length ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptMinLength"]),
          " -e ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptErrorRate"]),
          " --untrimmed-output ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.untrimmed.fastq")),
          " --untrimmed-paired-output ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt.untrimmed.fastq")),
          " -o ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1-{name}.fastq")),
          " -p ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt1-{name}.fastq")),
          " ",
          file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
          " ",
          file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[['exp_design']][i,"pair2"]),
          " > ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1.stdout")),
          " 2> ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1.stderr"))))
        #New read1 file
        temp_out = system(paste0(
          "cat ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt1-reverse.fastq")),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1-forward.fastq")),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.untrimmed.fastq")),
          " > ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt2"))))
        #New read2 file
        temp_out = system(paste0(
          "cat ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1-reverse.fastq")),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt1-forward.fastq")),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt.untrimmed.fastq")),
          " > ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt2"))))
        #Run cutadapt on the swapped  
        temp_out = system(paste0(
          "cutadapt",
          temp_options,
          " --discard-untrimmed --minimum-length ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptMinLength"]),
          " -e ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptErrorRate"]),
          " -j ",
          dimsum_meta[['numCores']],
          " -o ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt")),
          " -p ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt")),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt2")),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt2")),
          " > ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.stdout")),
          " 2> ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.stderr")))) 
      #Stranded library   
      }else{
        temp_out = system(paste0(
          "cutadapt",
          temp_options,
          temp_cut_options,
          " --discard-untrimmed --minimum-length ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptMinLength"]),
          " -e ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptErrorRate"]),
          " -j ",
          dimsum_meta[['numCores']],
          " -o ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt")),
          " -p ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt")),
          " ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
          " ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"], dimsum_meta[['exp_design']][i,"pair2"]),
          " > ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.stdout")),
          " 2> ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.stderr"))))
      }
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  temp_suffix <- ".cutadapt"
  dimsum_meta_new[["exp_design"]][,"pair1"] = paste0(dimsum_meta_new[["exp_design"]][,"pair1"], temp_suffix)
  dimsum_meta_new[["exp_design"]][,"pair2"] = paste0(dimsum_meta_new[["exp_design"]][,"pair2"], temp_suffix)
  dimsum_meta_new[['exp_design']][,"pair_directory"] <- cutadapt_outpath
  #Generate cutadapt report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_cutadapt_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
    return(dimsum_meta_new_report)
  }else{
    return(dimsum_meta_new)
  }
}

