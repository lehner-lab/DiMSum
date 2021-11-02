
#' dimsum__check_countfile
#'
#' Check whether user-specified count file correctly formatted.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#'
#' @return Reformatted data.table
#' @export
#' @import data.table
dimsum__check_countfile <- function(
  dimsum_meta,
  input_dt
  ){

  ### Nucleotide sequence checks (nt_seq column)
  #Check if mandatory columns present
  mandatory_cols <- c("nt_seq")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(input_dt)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from file specified by countPath ('nt_seq')"), call. = FALSE)
  }
  #Check nucleotide sequence column is of type character 
  if(sapply(input_dt, typeof)["nt_seq"]!="character"){
    stop("One or more invalid 'nt_seq' values in variant count file specified by 'countPath'. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }
  #Set nucleotide sequence to lower case
  input_dt[, nt_seq := tolower(nt_seq)]
  #Check nucleotide sequences are valid (ACGT characters only)
  if(sum(!input_dt[,unique(unlist(strsplit(nt_seq, "")))] %in% c('a', 'c', 'g', 't'))!=0){
    stop("One or more invalid 'nt_seq' values in variant count file specified by 'countPath'. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  ### Count column checks 
  #Check if sample name columns present
  mandatory_cols <- unique(dimsum_meta[["exp_design"]][,"sample_name"])
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(input_dt)))==FALSE)!=0){
    stop(paste0("One or more sample names in experimentDesign file missing from column names in variant count file specified by 'countPath'"), call. = FALSE)
  }
  #Check all count columns are of type integer 
  typeof_cols <- sapply(input_dt[,.SD,,.SDcols = names(input_dt)[names(input_dt)!="nt_seq"]], typeof)
  if(sum(typeof_cols!="integer")!=0){
    stop(paste0("Invalid type of sample count column in variant count file specified by 'countPath'. Only positive integers allowed (zero inclusive)."), call. = FALSE)
  }
  #Check all count columns positive integer zero inclusive
  if(input_dt[,min(.SD, na.rm = T),,.SDcols = names(input_dt)[names(input_dt)!="nt_seq"]]<0){
    stop(paste0("Invalid type of sample count column in variant count file specified by 'countPath'. Only positive integers allowed (zero inclusive)."), call. = FALSE)
  }

  ### Duplicated variants check
  if(input_dt[,sum(duplicated(nt_seq))]!=0){
    stop(paste0("Duplicated 'nt_seq' values not allowed in variant count file specified by 'countPath'."), call. = FALSE)
  }

  #Sample names (ignore 'technical_replicate' column)
  sample_names <- as.list(paste0(
    dimsum_meta[["exp_design"]][,"sample_name"], '_e', 
    dimsum_meta[["exp_design"]][,"experiment"], '_s', 
    dimsum_meta[["exp_design"]][,"selection_id"], '_b', 
    dimsum_meta[["exp_design"]][,"biological_replicate"], '_tNA_count', sep = ""))
  names(sample_names) <- dimsum_meta[["exp_design"]][,"sample_name"]

  #Reformat count column names
  names(input_dt)[names(input_dt)!="nt_seq"] <- unlist(sample_names[names(input_dt)[names(input_dt)!="nt_seq"]])

  return(input_dt)
}
