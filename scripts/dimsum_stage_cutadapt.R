
#dimsum_stage_cutadapt
#
# Run cutadapt on all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# cutadapt_outpath: cutadapt output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_cutadapt <- function(
  dimsum_meta,
  cutadapt_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL
  ){
  #Create/overwrite cutadapt directory (if executed)
  cutadapt_outpath <- gsub("/$", "", cutadapt_outpath)
  create_dimsum_dir(cutadapt_outpath, execute = execute, message = "DiMSum STAGE 5: CUTADAPT")  
  fastq_pair_list <- dimsum_meta[['exp_design']][,c('pair1', 'pair2')]
  rownames(fastq_pair_list) = 1:dim(fastq_pair_list)[1]
  #Cutadapt parameters specified?
  if( is.null(dimsum_meta[["cutadapt5First"]]) | is.null(dimsum_meta[["cutadapt5Second"]]) ){
    message("Skipping this stage (all cutadapt arguments need to be specified)")
    return(dimsum_meta)
  }else{
    #Options for removing constant regions from beginning or end of either read in pair
    temp_options = paste0(' -g ', dimsum_meta[["cutadapt5First"]], ' -G ', dimsum_meta[["cutadapt5Second"]])
    if( !is.null(dimsum_meta[["cutadapt3First"]]) ){temp_options = paste0(temp_options, " -a ", dimsum_meta[["cutadapt3First"]])}
    if( !is.null(dimsum_meta[["cutadapt3Second"]]) ){temp_options = paste0(temp_options, " -A ", dimsum_meta[["cutadapt3Second"]])}
    if( dimsum_meta[["cutadaptDiscardUntrimmed"]] ){temp_options = paste0(temp_options, " --discard-untrimmed ")}
    #Options for swapping read1 and read2
    temp_options_swap = paste0(' -g forward=', dimsum_meta[["cutadapt5First"]], ' -g reverse=', dimsum_meta[["cutadapt5Second"]])
    #Options for removing a fixed number of bases from beginning or end of either read in pair
    temp_cut_options = ''
    if( !is.null(dimsum_meta[["cutadaptCut5First"]]) ){temp_cut_options = paste0(temp_cut_options, " -u ", dimsum_meta[["cutadaptCut5First"]])}
    if( !is.null(dimsum_meta[["cutadaptCut3First"]]) ){temp_cut_options = paste0(temp_cut_options, " -u ", -dimsum_meta[["cutadaptCut3First"]])}
    if( !is.null(dimsum_meta[["cutadaptCut5Second"]]) ){temp_cut_options = paste0(temp_cut_options, " -U ", dimsum_meta[["cutadaptCut5Second"]])}
    if( !is.null(dimsum_meta[["cutadaptCut3Second"]]) ){temp_cut_options = paste0(temp_cut_options, " -U ", -dimsum_meta[["cutadaptCut3Second"]])}
    #Trim FASTQ file pairs
    message("Trimming FASTQ files with cutadapt:")
    all_fastq <- file.path(dimsum_meta[["exp_design"]]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
    print(all_fastq)
    message("Processing...")
    for(pair_name in rownames(fastq_pair_list)){
      #TODO: cutadapt binary path specifiable on commandline?
      #TEMP: don't trim files
      print(fastq_pair_list[pair_name,])
      #Check if this system command should be executed
      if(execute){
        #Not stranded library
        if( !dimsum_meta[["stranded"]] ){
          #Swap reads in fastq files according to adapter presence
          temp_out = system(paste0(
            "cutadapt",
            temp_options_swap,
            temp_cut_options,
            " --no-trim ",
            " --minimum-length ",
            as.character(dimsum_meta[["cutadaptMinLength"]]),
            " -e ",
            as.character(dimsum_meta[["cutadaptErrorRate"]]),
            " --untrimmed-output ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.untrimmed.fastq")),
            " --untrimmed-paired-output ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt.untrimmed.fastq")),
            " -o ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1-{name}.fastq")),
            " -p ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt1-{name}.fastq")),
            " ",
            file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][1]),
            " ",
            file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][2]),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1.stderr"))))
          #New read1 file
          temp_out = system(paste0(
            "cat ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt1-reverse.fastq")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1-forward.fastq")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.untrimmed.fastq")),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt2"))))
          #New read2 file
          temp_out = system(paste0(
            "cat ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1-reverse.fastq")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt1-forward.fastq")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt.untrimmed.fastq")),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt2"))))
          #Run cutadapt on the swapped  
          temp_out = system(paste0(
            "cutadapt",
            temp_options,
            " --minimum-length ",
            as.character(dimsum_meta[["cutadaptMinLength"]]),
            " -e ",
            as.character(dimsum_meta[["cutadaptErrorRate"]]),
            " -j ",
            num_cores,
            " -o ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt")),
            " -p ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt2")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt2")),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.stderr")))) 
        #Stranded library   
        }else{
          temp_out = system(paste0(
            "cutadapt",
            temp_options,
            temp_cut_options,
            " --minimum-length ",
            as.character(dimsum_meta[["cutadaptMinLength"]]),
            " -e ",
            as.character(dimsum_meta[["cutadaptErrorRate"]]),
            " -j ",
            num_cores,
            " -o ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt")),
            " -p ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt")),
            " ",
            file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][1]),
            " ",
            file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][2]),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.stderr"))))
        }
      }
    }
    #New experiment metadata
    dimsum_meta_new <- dimsum_meta
    #Update fastq metadata
    temp_suffix <- ".cutadapt"
    dimsum_meta_new[["exp_design"]]$pair1 = paste0(dimsum_meta_new[["exp_design"]]$pair1, temp_suffix)
    dimsum_meta_new[["exp_design"]]$pair2 = paste0(dimsum_meta_new[["exp_design"]]$pair2, temp_suffix)
    dimsum_meta_new[['exp_design']]$pair_directory <- cutadapt_outpath
    #Generate cutadapt report
    if(report){
      dimsum_meta_new_report <- dimsum_stage_cutadapt_report(dimsum_meta_new, report_outpath)
      return(dimsum_meta_new_report)
    }else{
      return(dimsum_meta_new)
    }
  }
}

