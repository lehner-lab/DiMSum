
#' dimsum__get_cutadapt_options
#'
#' Get command-line options for cutadapt.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param exp_design_row row number of exp_design table (required)
#' @param option_type one of default/cut/swap (default:"default")
#'
#' @return A string of command-line options for cutadapt
#' @export
dimsum__get_cutadapt_options <- function(
  dimsum_meta,
  exp_design_row,
  option_type = "default"){

  i <- exp_design_row

  #Options for removing constant regions from beginning or end of either read in pair
  if(option_type=="default"){
    temp_options <- ''
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt5First"]) ){temp_options <- paste0(' -g "', dimsum_meta[['exp_design']][i,"cutadapt5First"], '"')}
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt5Second"]) ){temp_options <- ifelse(dimsum_meta[['paired']], paste0(temp_options, ' -G "', dimsum_meta[['exp_design']][i,"cutadapt5Second"], '"'), temp_options)}
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt3First"]) ){temp_options <- paste0(temp_options, ' -a "', dimsum_meta[['exp_design']][i,"cutadapt3First"], '"')}
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt3Second"]) & dimsum_meta[['paired']] ){temp_options <- paste0(temp_options, ' -A "', dimsum_meta[['exp_design']][i,"cutadapt3Second"], '"')}
    #Discard untrimmed reads if at least one constant region trimming option specified
    if( !dimsum_meta[['exp_design']][i,"run_cutadapt_cutonly"] ){temp_options <- paste0(temp_options, " --discard-untrimmed")}
    return(temp_options)
  }

  #Options for removing a fixed number of bases from beginning or end of either read in pair
  if(option_type=="cut"){
    temp_cut_options <- ''
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadaptCut5First"]) ){temp_cut_options <- paste0(temp_cut_options, " -u ", dimsum_meta[['exp_design']][i,"cutadaptCut5First"])}
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadaptCut3First"]) ){temp_cut_options <- paste0(temp_cut_options, " -u ", -dimsum_meta[['exp_design']][i,"cutadaptCut3First"])}
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadaptCut5Second"]) & dimsum_meta[['paired']] ){temp_cut_options <- paste0(temp_cut_options, " -U ", dimsum_meta[['exp_design']][i,"cutadaptCut5Second"])}
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadaptCut3Second"]) & dimsum_meta[['paired']] ){temp_cut_options <- paste0(temp_cut_options, " -U ", -dimsum_meta[['exp_design']][i,"cutadaptCut3Second"])}
    return(temp_cut_options)
  }

  #Options for swapping read1 and read2
  if(option_type=="swap"){
    temp_options_swap <- ''
    if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt5First"]) & !is.na(dimsum_meta[['exp_design']][i,"cutadapt5Second"]) ){
      temp_options_swap <- paste0(' -g forward="', dimsum_meta[['exp_design']][i,"cutadapt5First"], '" -g reverse="', dimsum_meta[['exp_design']][i,"cutadapt5Second"], '"')
    }else if( !is.na(dimsum_meta[['exp_design']][i,"cutadapt3First"]) & !is.na(dimsum_meta[['exp_design']][i,"cutadapt3Second"]) ){
      temp_options_swap <- paste0(' -a forward="', dimsum_meta[['exp_design']][i,"cutadapt3First"], '" -a reverse="', dimsum_meta[['exp_design']][i,"cutadapt3Second"], '"')
    }
    return(temp_options_swap)
  }

}
