
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
  this_stage <- 2
  execute <- (dimsum_meta[["startStage"]] <= this_stage & dimsum_meta[["stopStage"]] >= this_stage)

  #WRAP not run or this stage after stopStage
  if(!is.null(dimsum_meta[["countPath"]]) | dimsum_meta[["stopStage"]] < this_stage){
    return(dimsum_meta)
  }

  #Save current workspace for debugging purposes
  if(save_workspace){dimsum__save_metadata(dimsum_meta = dimsum_meta, n = 2)}

  #Create/overwrite cutadapt directory (if executed)
  cutadapt_outpath <- gsub("/$", "", cutadapt_outpath)
  dimsum__create_dir(cutadapt_outpath, execute = execute, message = "DiMSum STAGE 2 (WRAP): TRIM CONSTANT REGIONS")  

  #If not trans library: convert to linked adapters if 3' constant region expected to be sequenced
  if(!dimsum_meta[["transLibrary"]]){
    dimsum_meta <- dimsum__convert_linked_adapters(dimsum_meta = dimsum_meta)
  }

  #Input files
  fastq_pair_list_nunique <- dimsum_meta[['exp_design']][,c('pair1', 'pair2')]
  all_fastq <- file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], unique(c(dimsum_meta[['exp_design']][,"pair1"], dimsum_meta[['exp_design']][,"pair2"])))
  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = all_fastq,
    stage_number = this_stage,
    execute = execute)

  #Not stranded library (and paired-end)
  if(!dimsum_meta[["stranded"]] & dimsum_meta[["paired"]]){
    #Swap reads in fastq files according to adapter presence
    dimsum__status_message("Swapping reads in FASTQ files according to adapter presence\n")
    if(execute){
      dimsum__swap_reads(dimsum_meta = dimsum_meta, cutadapt_outpath = cutadapt_outpath)
    }
  }

  #Check if all input files exist
  if(!dimsum_meta[["stranded"]] & dimsum_meta[["paired"]]){
    #Input files (after swapping)
    swap_fastq <- file.path(cutadapt_outpath, unique(c(paste0(dimsum_meta[['exp_design']][,"pair1"], ".cutadapt2.gz"), paste0(dimsum_meta[['exp_design']][,"pair1"], ".cutadapt2.gz"))))
    dimsum__check_files_exist(
      required_files = swap_fastq,
      stage_number = this_stage,
      execute = execute)
  }

  #Trim FASTQ file pairs
  dimsum__status_message("Trimming FASTQ files with cutadapt:\n")
  dimsum__status_message(paste0(all_fastq, "\n"))
  dimsum__status_message("Processing...\n")
  for(i in (1:dim(dimsum_meta[['exp_design']])[1])[!duplicated(fastq_pair_list_nunique)]){
    dimsum__status_message(paste0("\t", paste0(unique(unlist(dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')])), collapse = "\t"), "\n"))
    #Check if this system command should be executed
    if(execute){
      #Options for removing constant regions from beginning or end of either read in pair
      temp_options <- dimsum__get_cutadapt_options(dimsum_meta = dimsum_meta, exp_design_row = i)
      #Options for removing a fixed number of bases from beginning or end of either read in pair
      temp_cut_options <- dimsum__get_cutadapt_options(dimsum_meta = dimsum_meta, exp_design_row = i, option_type = "cut")
      #Insufficient cutadapt options specified to run cutadapt
      if(!dimsum_meta[['exp_design']][i,"run_cutadapt"]){
        #Copy files  
        file.copy(
          from = file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
          to = file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz")))
        if(dimsum_meta[["paired"]]){
          file.copy(
            from = file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], dimsum_meta[['exp_design']][i,"pair2"]),
            to = file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt.gz")))
          #Total number of FASTQ records
          temp_out <- system(paste0(
            "wc -l ",
            file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
            " ",
            file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], dimsum_meta[['exp_design']][i,"pair2"]),
            " > ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz.stdout")))) 
        }else{
          #Total number of FASTQ records
          temp_out <- system(paste0(
            "wc -l ",
            file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
            " > ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz.stdout")))) 
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
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz")),
          " -p ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt.gz")),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt2.gz")),
          " ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt2.gz")),
          " > ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz.stdout")),
          " 2> ",
          file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz.stderr")))) 
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
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz")),
            " -p ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt.gz")),
            " ",
            file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
            " ",
            file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], dimsum_meta[['exp_design']][i,"pair2"]),
            " > ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz.stderr"))))
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
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz")),
            " ",
            file.path(dimsum_meta[["exp_design"]][1,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
            " > ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.gz.stderr"))))          
        }
      }
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  temp_suffix <- ".cutadapt.gz"
  dimsum_meta_new[["exp_design"]][,"pair1"] <- paste0(dimsum_meta_new[["exp_design"]][,"pair1"], temp_suffix)
  dimsum_meta_new[["exp_design"]][,"pair2"] <- paste0(dimsum_meta_new[["exp_design"]][,"pair2"], temp_suffix)
  dimsum_meta_new[['exp_design']][,"pair_directory"] <- cutadapt_outpath
  #Delete files when last stage complete
  if(!dimsum_meta_new[["retainIntermediateFiles"]]){
    if(dimsum_meta_new[["stopStage"]]==this_stage){
      if(!is.null(dimsum_meta_new[["deleteIntermediateFiles"]])){suppressWarnings(temp_out <- file.remove(dimsum_meta_new[["deleteIntermediateFiles"]]))}
      if(!is.null(dimsum_meta_new[["touchIntermediateFiles"]])){suppressWarnings(temp_out <- file.create(dimsum_meta_new[["touchIntermediateFiles"]]))}
    }else{
      dimsum_meta_new[["deleteIntermediateFiles"]] <- c(dimsum_meta_new[["deleteIntermediateFiles"]], 
        file.path(cutadapt_outpath, dir(cutadapt_outpath, "*.fastq$")),
        file.path(cutadapt_outpath, dir(cutadapt_outpath, "*.cutadapt.gz$")),
        file.path(cutadapt_outpath, dir(cutadapt_outpath, "*.cutadapt2.gz$")))
    }
  }
  #Generate cutadapt report
  if(report){
    tryCatch({
      dimsum_meta_new_report <- dimsum__cutadapt_report(dimsum_meta = dimsum_meta_new, report_outpath = report_outpath)
      }, error=function(e){
        dimsum__status_message("There were problems while running 'dimsum__cutadapt_report'\n")
        dimsum_meta_new_report <- dimsum_meta_new
        })
    return(dimsum_meta_new_report)
  }else{
    return(dimsum_meta_new)
  }
}

