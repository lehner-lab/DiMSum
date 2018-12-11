
#get_experiment_metadata
#
# Get and format metadata from experiment design file.
#
# experiment_design_path: path to experiment design file (required)
#
# Returns: a data frame
#
get_experiment_design <- function(
  dimsum_meta){
  #Read experimental design
  if(!file.exists(dimsum_meta[["experiment_design_path"]])){
    stop("Experiment design file not found. Please check that the --experimentDesignPath (-e) command-line option is correctly set.", call. = FALSE)
  }
  exp_design <- read.table(dimsum_meta[["experiment_design_path"]], header = T, stringsAsFactors = F, sep="\t")
  #Add original FASTQ directory
  if(!file.exists(dimsum_meta[["fastq_path_original"]])){
    stop("Input FASTQ directory not found. Please check that the --fastqFileDir (-l) command-line option is correctly set.", call. = FALSE)
  }
  exp_design[,"pair_directory"] <- dimsum_meta[["fastq_path_original"]]
  #Add sample specific cutadapt options
  if((!"cutadapt5First" %in% colnames(exp_design) & is.null(dimsum_meta[["cutadapt5First"]])) | (!"cutadapt5Second" %in% colnames(exp_design) & is.null(dimsum_meta[["cutadapt5Second"]]))){
    stop("Sequence of 5' constant region not found. Please check that the --cutadapt5First and --cutadapt5Second command-line options are correctly set.", call. = FALSE)
  }
  if(!"cutadaptCut5First" %in% colnames(exp_design)){exp_design[,"cutadaptCut5First"] <- ifelse(is.null(dimsum_meta[["cutadaptCut5First"]]), NA, dimsum_meta[["cutadaptCut5First"]])}
  if(!"cutadaptCut5Second" %in% colnames(exp_design)){exp_design[,"cutadaptCut5Second"] <- ifelse(is.null(dimsum_meta[["cutadaptCut5Second"]]), NA, dimsum_meta[["cutadaptCut5Second"]])}
  if(!"cutadaptCut3First" %in% colnames(exp_design)){exp_design[,"cutadaptCut3First"] <- ifelse(is.null(dimsum_meta[["cutadaptCut3First"]]), NA, dimsum_meta[["cutadaptCut3First"]])}
  if(!"cutadaptCut3Second" %in% colnames(exp_design)){exp_design[,"cutadaptCut3Second"] <- ifelse(is.null(dimsum_meta[["cutadaptCut3Second"]]), NA, dimsum_meta[["cutadaptCut3Second"]])}
  if(!"cutadapt5First" %in% colnames(exp_design)){exp_design[,"cutadapt5First"] <- ifelse(is.null(dimsum_meta[["cutadapt5First"]]), NA, dimsum_meta[["cutadapt5First"]])}
  if(!"cutadapt5Second" %in% colnames(exp_design)){exp_design[,"cutadapt5Second"] <- ifelse(is.null(dimsum_meta[["cutadapt5Second"]]), NA, dimsum_meta[["cutadapt5Second"]])}
  if(!"cutadapt3First" %in% colnames(exp_design)){exp_design[,"cutadapt3First"] <- ifelse(is.null(dimsum_meta[["cutadapt3First"]]), NA, dimsum_meta[["cutadapt3First"]])}
  if(!"cutadapt3Second" %in% colnames(exp_design)){exp_design[,"cutadapt3Second"] <- ifelse(is.null(dimsum_meta[["cutadapt3Second"]]), NA, dimsum_meta[["cutadapt3Second"]])}
  if(!"cutadaptMinLength" %in% colnames(exp_design)){exp_design[,"cutadaptMinLength"] <- ifelse(is.null(dimsum_meta[["cutadaptMinLength"]]), NA, dimsum_meta[["cutadaptMinLength"]])}
  if(!"cutadaptErrorRate" %in% colnames(exp_design)){exp_design[,"cutadaptErrorRate"] <- ifelse(is.null(dimsum_meta[["cutadaptErrorRate"]]), NA, dimsum_meta[["cutadaptErrorRate"]])}
  if(!"cutadaptDiscardUntrimmed" %in% colnames(exp_design)){exp_design[,"cutadaptDiscardUntrimmed"] <- ifelse(is.null(dimsum_meta[["cutadaptDiscardUntrimmed"]]), NA, dimsum_meta[["cutadaptDiscardUntrimmed"]])}
  return(exp_design)
}
