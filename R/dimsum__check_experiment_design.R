
#' dimsum__check_experiment_design
#'
#' Validate metadata from experiment design file.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param exp_design experiment design data.frame (required)
#'
#' @return Nothing
#' @export
dimsum__check_experiment_design <- function(
  dimsum_meta,
  exp_design
  ){

  ### Column existence checks 
  #Check if mandatory columns present
  mandatory_cols <- c("sample_name", "experiment_replicate", "selection_id", "selection_replicate", "technical_replicate", "pair1", "pair2")
  if(sum(unlist(lapply(mandatory_cols, "%in%", colnames(exp_design)))==FALSE)!=0){
    stop(paste0("One or more mandatory columns missing from experimentDesign file ('sample_name', 'experiment_replicate', 'selection_id', 'selection_replicate', 'technical_replicate', 'pair1', 'pair2')"), call. = FALSE)
  }

  ### FASTQ file checks (pair1/pair2 columns)
  #Check file name columns are of type character 
  if(typeof(exp_design[,"pair1"])!="character" | typeof(exp_design[,"pair2"])!="character"){
    stop(paste0("One or more invalid FASTQ file name values in experimentDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for duplicated FASTQ files
  if(sum(duplicated(exp_design[,c("pair1", "pair2")]))!=0 & !dimsum_meta[["experimentDesignPairDuplicates"]]){
    stop(paste0("Duplicate FASTQ files not allowed in experimentDesign file columns: 'pair1' and 'pair2'"), call. = FALSE)
  }

  ### Sample name checks (sample_name column)
  #Check sample name column is of type character 
  if(typeof(exp_design[,"sample_name"])!="character"){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (only characters allowed)"), call. = FALSE)
  }
  #Check for duplicated sample names (after collapsing on technical replicate id)
  sample_names_collapsed <- unique(exp_design[,c("sample_name", "experiment_replicate", "selection_id", "selection_replicate")])[,"sample_name"]
  if(sum(duplicated(sample_names_collapsed))!=0){
    stop(paste0("Duplicate 'sample_name' values not allowed in experimentDesign file (if not technical replicates)"), call. = FALSE)
  }
  #Check only alphanumeric characters
  if(sum(!grepl("^[A-Za-z0-9]+$", exp_design[,"sample_name"], perl = T))!=0){
    stop(paste0("One or more invalid 'sample_name' values in experimentDesign file (only alphanumeric characters allowed)"), call. = FALSE)
  }
  #Check that all input sample names for the same experiment_replicate are identical
  experiment_input <- unique(exp_design[exp_design[,"selection_id"]==0,c("sample_name", "experiment_replicate")])
  if(sum(duplicated(experiment_input[,"experiment_replicate"]))!=0){
    stop(paste0("One or more invalid input 'sample_name' values in experimentDesign file (must be identical for the same experiment replicate)"), call. = FALSE)
  }
  #Check that all output sample names for the same experiment_replicate and selection_replicate are identical
  experiment_output <- unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "experiment_replicate", "selection_replicate")])
  if(sum(duplicated(experiment_output[,c("experiment_replicate", "selection_replicate")]))!=0){
    stop(paste0("One or more invalid output 'sample_name' values in experimentDesign file (must be identical for the same experiment replicate and selection replicate)"), call. = FALSE)
  }

  ### Experiment id checks (experiment_replicate column)
  #Check experiment_replicate strictly positive integer
  if(typeof(exp_design[,"experiment_replicate"])!="integer"){
    stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"experiment_replicate"]<=0)!=0){
    stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  # #Check that at least two experimental replicates exist (requirement to run error model)
  # if(length(unique(exp_design[,"experiment_replicate"]))<2){
  #   stop(paste0("One or more invalid 'experiment_replicate' values in experimentDesign file (at least two replicates required)"), call. = FALSE)
  # }

  ### Selection id checks (selection_id column)
  #Check selection_id strictly positive integer (zero inclusive)
  if(typeof(exp_design[,"selection_id"])!="integer"){
    stop(paste0("One or more invalid 'selection_id' values in experimentDesign file (only 0 or positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"selection_id"]<0)!=0){
    stop(paste0("One or more invalid 'selection_id' values in experimentDesign file (only 0 or positive integers allowed)"), call. = FALSE)
  }

  ### Minimum required matrix data checks
  #Check if at least one input and output sample present for each experiment_replicate
  experiment_input <- unique(exp_design[exp_design[,"selection_id"]==0,"experiment_replicate"])
  experiment_output <- unique(exp_design[exp_design[,"selection_id"]!=0,"experiment_replicate"])
  if(length(experiment_input)!=length(experiment_output) | sum(!experiment_input %in% experiment_output)!=0){
    stop(paste0("All experiment replicates should have at least one sample before and after selection in experimentDesign file"), call. = FALSE)
  }

  ### Biological replicate id checks (selection_replicate column)
  #Check that selection_replicate column blank (NA) for input samples
  if(sum(!is.na(exp_design[exp_design[,"selection_id"]==0,"selection_replicate"]))!=0){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (leave blank for input samples)"), call. = FALSE)
  }
  #Check selection_replicate strictly positive integer
  if(typeof(exp_design[,"selection_replicate"])!="integer"){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
  }
  if(sum(exp_design[,"selection_replicate"]<=0, na.rm = T)!=0){
    stop(paste0("One or more invalid 'selection_replicate' values in experimentDesign file (only positive integers allowed)"), call. = FALSE)
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
  #Check for duplicated rows in the following columns: "experiment_replicate", "selection_id", "selection_replicate", "technical_replicate"
  if(sum(duplicated(exp_design[,c("experiment_replicate", "selection_id", "selection_replicate", "technical_replicate")]))!=0){
    stop(paste0("One or more duplicated rows in experimentDesign file matrix (sample rows should be unique)"), call. = FALSE)
  }

  # ### 5'/3' adapter/cut existence checks
  # #Check that at least one cutadapt argument specified for all samples
  # if(sum(is.na(exp_design[,"cutadapt5First"]))!=0 & sum(is.na(exp_design[,"cutadapt5Second"]))!=0 & sum(is.na(exp_design[,"cutadapt3First"]))!=0 & sum(is.na(exp_design[,"cutadapt3Second"]))!=0){
  #   if(sum(is.na(exp_design[,"cutadaptCut5First"]))!=0 & sum(is.na(exp_design[,"cutadaptCut5Second"]))!=0 & sum(is.na(exp_design[,"cutadaptCut3First"]))!=0 & sum(is.na(exp_design[,"cutadaptCut3Second"]))!=0){
  #     stop("At least one of cutadapt5First, cutadapt5Second, cutadapt3First, cutadapt3Second, cutadaptCut5First, cutadaptCut5Second, cutadaptCut3First, cutadaptCut3Second must be specified for all samples. Please check that the corresponding experiment design file columns are correct.", call. = FALSE)
  #   }
  # }

  ### Constant region checks
  #Check constant region columns are of type character (or logical i.e. all empty/NA)
  if(sum(!unlist(lapply(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")], typeof)) %in% c("character", "logical"))!=0){
    stop(paste0("One or more invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G)."), call. = FALSE)
  }
  #Check constant region sequences are valid (ACGT characters only)
  all_constant_regions <- as.character(unlist(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")]))
  all_constant_regions <- gsub(";required\\.\\.\\.|;optional\\.\\.\\.", "\\.\\.\\.", all_constant_regions)
  all_constant_regions <- gsub(";required$|;optional$", "", all_constant_regions)
  all_characters <- unique(unlist(strsplit(gsub("[ACGT]\\.\\.\\.[ACGT]", "", all_constant_regions), "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T", NA))!=0){
    stop("Invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }

  ### Misc cutadapt argument checks
  #Check cutadaptCut... columns are of type integer (or logical i.e. all empty/NA)
  if(sum(!unlist(lapply(exp_design[,c("cutadaptCut5First", "cutadaptCut5Second", "cutadaptCut3First", "cutadaptCut3Second")], typeof)) %in% c("integer", "logical"))!=0){
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
  #Check cutadaptErrorRate double
  if(typeof(exp_design[,"cutadaptErrorRate"])!="double"){
    stop("Invalid 'cutadaptErrorRate' argument. Only positive doubles less than 1 allowed (zero inclusive).", call. = FALSE)
  }
  #Check cutadaptErrorRate positive double less than 1
  if(sum(exp_design[,c("cutadaptErrorRate")]<0)!=0 | sum(exp_design[,c("cutadaptErrorRate")]>=1)!=0){
    stop("Invalid 'cutadaptErrorRate' argument. Only positive doubles less than 1 allowed (zero inclusive).", call. = FALSE)
  }
  #Check cutadaptOverlap integer
  if(typeof(exp_design[,"cutadaptOverlap"])!="integer"){
    stop("Invalid 'cutadaptOverlap' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }
  #Check cutadaptOverlap positive integer zero inclusive
  if(sum(exp_design[,c("cutadaptOverlap")]<0)!=0){
    stop("Invalid 'cutadaptOverlap' argument. Only positive integers allowed (zero inclusive).", call. = FALSE)
  }

  ### Generations checks
  #Check generations column is of type double (or logical i.e. all empty/NA)
  if(!typeof(exp_design[,"generations"]) %in% c("double", "logical")){
    stop(paste0("One or more invalid generations values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check generations column is strictly positive
  if(sum(exp_design[,"generations"]<=0, na.rm=T)!=0){
    stop(paste0("One or more invalid generations values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check that all output samples have a generations value (if one or more specified)
  if(typeof(exp_design[,"generations"])=="double" & sum(is.na(exp_design[exp_design[,"selection_id"]==1,"generations"]))!=0){
    stop(paste0("One or more missing generations values. Generations values must be specified for all output samples (or none)."), call. = FALSE)
  }
  #Check that generations values are identical for technical output replicates (if one or more specified)
  temp_n_output_replicates <- length(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name")]))
  if(typeof(exp_design[,"generations"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "generations")]))[1]!=temp_n_output_replicates){
    stop(paste0("Generations values not identical for technical output replicates. Generations values must be specified for all output samples (or none)."), call. = FALSE)
  }

  ### Selection time checks
  #Check selection_time column is of type double (or logical i.e. all empty/NA)
  if(!typeof(exp_design[,"selection_time"]) %in% c("double", "logical")){
    stop(paste0("One or more invalid selection_time values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check selection_time column is strictly positive
  if(sum(exp_design[,"selection_time"]<=0, na.rm=T)!=0){
    stop(paste0("One or more invalid selection_time values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check that all output samples have a selection_time value (if one or more specified)
  if(typeof(exp_design[,"selection_time"])=="double" & sum(is.na(exp_design[exp_design[,"selection_id"]==1,"selection_time"]))!=0){
    stop(paste0("One or more missing selection_time values. Selection time values must be specified for all output samples (or none)."), call. = FALSE)
  }
  #Check that selection_time values are identical for technical output replicates (if one or more specified)
  temp_n_output_replicates <- length(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name")]))
  if(typeof(exp_design[,"selection_time"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "selection_time")]))[1]!=temp_n_output_replicates){
    stop(paste0("Selection time values not identical for technical output replicates. Selection time values must be specified for all output samples (or none)."), call. = FALSE)
  }

  ### Cell density checks
  #Check cell_density column is of type double (or logical i.e. all empty/NA)
  if(!typeof(exp_design[,"cell_density"]) %in% c("double", "logical")){
    stop(paste0("One or more invalid cell_density values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check cell_density column is strictly positive
  if(sum(exp_design[,"cell_density"]<=0, na.rm=T)!=0){
    stop(paste0("One or more invalid cell_density values. Only positive doubles allowed (zero exclusive)."), call. = FALSE)
  }
  #Check that all samples have a cell_density value (if one or more specified)
  if(typeof(exp_design[,"cell_density"])=="double" & sum(is.na(exp_design[,"cell_density"]))!=0){
    stop(paste0("One or more missing cell_density values. Cell density values must be specified for all samples (or none)."), call. = FALSE)
  }
  #Check that cell_density values are identical for all input technical replicates (if one or more specified)
  temp_n_input_replicates <- length(unique(exp_design[exp_design[,"selection_id"]==0,c("sample_name")]))
  if(typeof(exp_design[,"cell_density"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]==0,c("sample_name", "cell_density")]))[1]!=temp_n_input_replicates){
    stop(paste0("Cell density values not identical for technical input replicates. Cell density values must be specified for all samples (or none)."), call. = FALSE)
  }
  #Check that cell_density values are identical for all output technical replicates (if one or more specified)
  temp_n_output_replicates <- length(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name")]))
  if(typeof(exp_design[,"cell_density"])=="double" & dim(unique(exp_design[exp_design[,"selection_id"]!=0,c("sample_name", "cell_density")]))[1]!=temp_n_output_replicates){
    stop(paste0("Cell density values not identical for technical output replicates. Cell density values must be specified for all samples (or none)."), call. = FALSE)
  }

}


