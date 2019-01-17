
#' dimsum__check_experiment_design
#'
#' Validate metadata from experiment design file.
#'
#' @param exp_design experiment design data.frame (required)
#'
#' @return Nothing
#' @export
dimsum__check_experiment_design <- function(
  exp_design
  ){

  ### Column existence checks 
  #Check if mandatory columns present
  mandatory_cols <- c("sample_name", "experiment", "selection_id", "biological_replicate", "technical_replicate", "pair1", "pair2")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(exp_design)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from experimentDesign file ('sample_name', 'experiment', 'selection_id', 'biological_replicate', 'technical_replicate', 'pair1', 'pair2')"), call. = FALSE)
  }

  ### FASTQ file checks (pair1/pair2 columns)
  #Check file name columns are of type character 
  if(typeof(exp_design[,"pair1"])!="character" | typeof(exp_design[,"pair2"])!="character"){
    stop(paste0("One or more invalid FASTQ file name values in experimentDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for duplicated FASTQ files
  if(sum(duplicated(exp_design[,c("pair1", "pair2")]))!=0){
    stop(paste0("Duplicate FASTQ files not allowed in experimentDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }

  ### Sample name checks (sample_name column)
  #Check sample name column is of type character 
  if(typeof(exp_design[,"sample_name"])!="character"){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for duplicated sample names (after collapsing on technical replicate id)
  sample_names_collapsed <- unique(exp_design[,c("sample_name", "experiment", "selection_id", "biological_replicate")])[,"sample_name"]
  if(sum(duplicated(sample_names_collapsed))!=0){
    stop(paste0("Duplicate 'sample_name' values not allowed in experimentDesign file (if not technical replicates)"), call. = FALSE)
  }
  #Check only alphanumeric characters
  if(sum(!grepl("^[A-Za-z0-9]+$", exp_design[,"sample_name"], perl = T))!=0){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (only alphanumeric characters allowed)"), call. = FALSE)
  }

  ### Experiment id checks (experiment column)
  #Check experiment strictly positive integer
  if(typeof(exp_design[,"experiment"])!="integer"){
    stop(paste0("One or more invalid 'experiment' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"experiment"]<=0)!=0){
    stop(paste0("One or more invalid 'experiment' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }

  ### Selection id checks (selection_id column)
  #Check selection_id strictly positive integer (zero inclusive)
  if(typeof(exp_design[,"selection_id"])!="integer"){
    stop(paste0("One or more invalid 'selection_id' values in experimentDesign file (only 0 or positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"selection_id"]<0)!=0){
    stop(paste0("One or more invalid 'selection_id' values in experimentDesign file (only 0 or positive integers allowed)"), call. = FALSE)
  }

  ### Minimum required matrix data checks
  #Check if at least one input and output sample present for each experiment
  experiment_input <- unique(exp_design[exp_design[,"selection_id"]==0,"experiment"])
  experiment_output <- unique(exp_design[exp_design[,"selection_id"]!=0,"experiment"])
  if(length(experiment_input)!=length(experiment_output) | sum(!experiment_input %in% experiment_output)!=0){
    stop(paste0("All experiments should have at least one sample before and after selection in experimentDesign file"), call. = FALSE)
  }

  ### Biological replicate id checks (biological_replicate column)
  #Check that biological_replicate column blank (NA) for input samples
  if(sum(!is.na(exp_design[exp_design[,"selection_id"]==0,"biological_replicate"]))!=0){
    stop(paste0("One or more invalid 'biological_replicate' values in experimentDesign file (leave blank for input samples)"), call. = FALSE)
  }
  #Check biological_replicate strictly positive integer
  if(typeof(exp_design[,"biological_replicate"])!="integer"){
    stop(paste0("One or more invalid 'biological_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"biological_replicate"]<=0, na.rm = T)!=0){
    stop(paste0("One or more invalid 'biological_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }

  ### Technical replicate id checks (technical_replicate column)
  #Check technical_replicate strictly positive integer (or all missing/NA)
  if(typeof(exp_design[,"technical_replicate"])!="integer" & sum(!is.na(exp_design[,"technical_replicate"]))!=0){
    stop(paste0("One or more invalid 'technical_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"technical_replicate"]<=0, na.rm = T)!=0){
    stop(paste0("One or more invalid 'technical_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }

  ### Duplicate matrix row checks
  #Check for duplicated rows in the following columns: "experiment", "selection_id", "biological_replicate", "technical_replicate"
  if(sum(duplicated(exp_design[,c("experiment", "selection_id", "biological_replicate", "technical_replicate")]))!=0){
    stop(paste0("One or more duplicated rows in experimentDesign file matrix (sample rows should be unique)"), call. = FALSE)
  }

  ### 5' adapter existence checks
  #Check that each sample has a 5' adapter (constant region) specified
  if(sum(is.na(exp_design[,"cutadapt5First"]))!=0 | sum(is.na(exp_design[,"cutadapt5Second"]))!=0){
    stop("Sequence of 5' constant region not found for some samples. Please check that the corresponding experiment design file columns are correct.", call. = FALSE)
  }

  ### Constant region checks
  #Check constant region columns are of type character (or logical i.e. all empty/NA)
  if(sum(!unlist(lapply(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")], typeof)) %in% c("character", "logical"))!=0){
    stop(paste0("One or more invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G)."), call. = FALSE)
  }
  #Check constant region sequences are valid (ACGT characters only)
  all_characters <- unique(unlist(strsplit(unlist(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")]), "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T", NA))!=0){
    stop("Invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  ### Misc cutadapt argument checks
  #Check cutadaptCut... columns are of type integer (or logical i.e. all empty/NA)
  if(sum(!unlist(lapply(exp_design[,c("cutadaptCut5First", "cutadaptCut5Second", "cutadaptCut3First", "cutadaptCut3Second")], typeof)) %in% c("character", "logical"))!=0){
    stop(paste0("One or more invalid 'cutadaptCut...' arguments. Only positive integers allowed (zero exclusive)."), call. = FALSE)
  }
  #Check strictly positive integer cutadaptCut... arguments (if not NA)
  if(sum(unlist(exp_design[,c("cutadaptCut5First", "cutadaptCut5Second", "cutadaptCut3First", "cutadaptCut3Second")])<=0, na.rm = T)!=0){
    stop("Invalid 'cutadaptCut...' arguments. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }
  #Check cutadaptMinLength strictly positive integer
  if(typeof(exp_design[,"cutadaptMinLength"])!="integer"){
    stop("Invalid 'cutadaptMinLength' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }
  #Check cutadaptMinLength argument
  if(sum(exp_design[,c("cutadaptMinLength")]<1)!=0){
    stop("Invalid 'cutadaptMinLength' argument. Only positive integers allowed (zero exclusive).", call. = FALSE)
  }
  #Check cutadaptErrorRate strictly positive integer
  if(typeof(exp_design[,"cutadaptErrorRate"])!="double"){
    stop("Invalid 'cutadaptErrorRate' argument. Only positive doubles allowed (zero exclusive).", call. = FALSE)
  }
  #Check cutadaptErrorRate argument
  if(sum(exp_design[,c("cutadaptErrorRate")]<0)!=0){
    stop("Invalid 'cutadaptErrorRate' argument. Only positive doubles allowed (zero inclusive).", call. = FALSE)
  }

}


