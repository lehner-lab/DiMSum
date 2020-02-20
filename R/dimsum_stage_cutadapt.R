
#' dimsum_stage_cutadapt
#'
#' Run cutadapt on all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param cutadapt_outpath cutadapt output path (required)
#' @param report whether or not to generate cutadapt summary plots (default: TRUE)
#' @param report_outpath cutadapt report output path
#' @param save_workspace whether or not to save the current workspace (default: TRUE)
#'
#' @return an updated experiment metadata object
#' @export
dimsum_stage_cutadapt <- function(
  dimsum_meta,
  cutadapt_outpath,
  report = TRUE,
  report_outpath = NULL,
  save_workspace = TRUE
  ){
  #Whether or not to execute the system command
  this_stage <- 3
  execute <- (dimsum_meta[["startStage"]] <= this_stage & (dimsum_meta[["stopStage"]] == 0 | dimsum_meta[["stopStage"]] >= this_stage))
  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}
  #Create/overwrite cutadapt directory (if executed)
  cutadapt_outpath <- gsub("/$", "", cutadapt_outpath)
  dimsum__create_dir(cutadapt_outpath, execute = execute, message = "DiMSum STAGE 3: TRIM CONSTANT REGIONS")  
  fastq_pair_list_nunique <- dimsum_meta[['exp_design']][,c('pair1', 'pair2')]
  #If not trans library: convert to linked adapters if 3' constant region expected to be sequenced
  if(!dimsum_meta[["transLibrary"]]){
    dimsum_meta <- dimsum__convert_linked_adapters(dimsum_meta = dimsum_meta)
  }
  #Not stranded library (and paired-end)
  if(!dimsum_meta[["stranded"]] & dimsum_meta[["paired"]]){
    #Swap reads in fastq files according to adapter presence
    message("Swapping reads in FASTQ files according to adapter presence")
    if(execute){
      dimsum__swap_reads(dimsum_meta = dimsum_meta, cutadapt_outpath = cutadapt_outpath)
    }
  }
  #Trim FASTQ file pairs
  message("Trimming FASTQ files with cutadapt:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], unique(c(dimsum_meta[['exp_design']][,"pair1"], dimsum_meta[['exp_design']][,"pair2"])))
  print(all_fastq)
  message("Processing...")
  for(i in (1:dim(dimsum_meta[['exp_design']])[1])[!duplicated(fastq_pair_list_nunique)]){
    message(paste0("\t", unique(unlist(dimsum_meta[['exp_design']][i,c('pair1', 'pair2')]))))
    #Check if this system command should be executed
    if(execute){
      #Options for removing constant regions from beginning or end of either read in pair
      temp_options <- dimsum__get_cutadapt_options(dimsum_meta = dimsum_meta, exp_design_row = i)
      #Options for removing a fixed number of bases from beginning or end of either read in pair
      temp_cut_options <- dimsum__get_cutadapt_options(dimsum_meta = dimsum_meta, exp_design_row = i, option_type = "cut")
      #Insufficient cutadapt options specified to run cutadapt
      if(!dimsum_meta[['exp_design']][i,"run_cutadapt"]){
        #Copy files  
        temp_out <- system(paste0(
          "cp ",
          file.path(dimsum_meta[["exp_design"]][,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt"))))
        if(dimsum_meta[["paired"]]){
          temp_out <- system(paste0(
            "cp ",
            file.path(dimsum_meta[["exp_design"]][,"pair_directory"], dimsum_meta[['exp_design']][i,"pair2"]),
            " ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt"))))
          #Total number of FASTQ records
          temp_out <- system(paste0(
            "wc -l ",
            file.path(dimsum_meta[["exp_design"]][,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
            " ",
            file.path(dimsum_meta[["exp_design"]][,"pair_directory"], dimsum_meta[['exp_design']][i,"pair2"]),
            " > ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.stdout")))) 
        }else{
          #Total number of FASTQ records
          temp_out <- system(paste0(
            "wc -l ",
            file.path(dimsum_meta[["exp_design"]][,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
            " > ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.stdout")))) 
        }
      #Not stranded library (and paired-end)
      }else if( !dimsum_meta[["stranded"]] & dimsum_meta[["paired"]]){
        #Run cutadapt on the swapped  
        temp_out <- system(paste0(
          "cutadapt",
          temp_options,
          " --minimum-length ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptMinLength"]),
          " -e ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptErrorRate"]),
          " -O ",
          as.character(dimsum_meta[['exp_design']][i,"cutadaptOverlap"]),
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
        if(dimsum_meta[["paired"]]){
          temp_out <- system(paste0(
            "cutadapt",
            temp_options,
            temp_cut_options,
            " --minimum-length ",
            as.character(dimsum_meta[['exp_design']][i,"cutadaptMinLength"]),
            " -e ",
            as.character(dimsum_meta[['exp_design']][i,"cutadaptErrorRate"]),
            " -O ",
            as.character(dimsum_meta[['exp_design']][i,"cutadaptOverlap"]),
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
        }else{
          temp_out <- system(paste0(
            "cutadapt",
            temp_options,
            temp_cut_options,
            " --minimum-length ",
            as.character(dimsum_meta[['exp_design']][i,"cutadaptMinLength"]),
            " -e ",
            as.character(dimsum_meta[['exp_design']][i,"cutadaptErrorRate"]),
            " -O ",
            as.character(dimsum_meta[['exp_design']][i,"cutadaptOverlap"]),
            " -j ",
            dimsum_meta[['numCores']],
            " -o ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt")),
            " ",
            file.path(dimsum_meta[["exp_design"]][,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
            " > ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.stderr"))))          
        }
      }
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  temp_suffix <- ".cutadapt"
  dimsum_meta_new[["exp_design"]][,"pair1"] <- paste0(dimsum_meta_new[["exp_design"]][,"pair1"], temp_suffix)
  dimsum_meta_new[["exp_design"]][,"pair2"] <- paste0(dimsum_meta_new[["exp_design"]][,"pair2"], temp_suffix)
  dimsum_meta_new[['exp_design']][,"pair_directory"] <- cutadapt_outpath
  #Delete files when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      temp_out <- mapply(system, dimsum_meta_new[["deleteIntermediateFiles"]], MoreArgs = list(ignore.stdout = T, ignore.stderr = T))
    }else{
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], 
        paste0("rm ", file.path(cutadapt_outpath, "*.fastq")),
        paste0("rm ", file.path(cutadapt_outpath, "*.cutadapt")),
        paste0("rm ", file.path(cutadapt_outpath, "*.cutadapt2")))
    }
  }
  #Generate cutadapt report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_cutadapt_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
    return(dimsum_meta_new_report)
  }else{
    return(dimsum_meta_new)
  }
}

