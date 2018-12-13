
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
  #Convert empty string constant region sequences to NA
  exp_design[which(exp_design[,"cutadapt5First"]==""),"cutadapt5First"] <- NA
  exp_design[which(exp_design[,"cutadapt5Second"]==""),"cutadapt5Second"] <- NA
  exp_design[which(exp_design[,"cutadapt3First"]==""),"cutadapt3First"] <- NA
  exp_design[which(exp_design[,"cutadapt3Second"]==""),"cutadapt3Second"] <- NA
  #Check that each sample has a 5' adapter (constant region) specified
  if(sum(is.na(exp_design[,"cutadapt5First"]))!=0 | sum(is.na(exp_design[,"cutadapt5Second"]))!=0){
    stop("Sequence of 5' constant region not found for some samples. Please check that the corresponding experiment design file columns are correct.", call. = FALSE)
  }
  #Check constant region sequences are valid (ACGT characters only)
  all_characters <- unique(unlist(strsplit(unlist(exp_design[,c("cutadapt5First", "cutadapt5Second", "cutadapt3First", "cutadapt3Second")]), "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T", NA))!=0){
    stop("Invalid constant region sequences. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }
  #If not trans library: reverse complement cutadapt 5' constant regions to obtain 3' constant regions (if not already supplied)
  if(!dimsum_meta[["transLibrary"]]){
    exp_design[is.na(exp_design[,"cutadapt3First"]),"cutadapt3First"] <- as.character(reverseComplement(DNAStringSet(exp_design[is.na(exp_design[,"cutadapt3First"]),"cutadapt5Second"])))
    exp_design[is.na(exp_design[,"cutadapt3Second"]),"cutadapt3Second"] <- as.character(reverseComplement(DNAStringSet(exp_design[is.na(exp_design[,"cutadapt3Second"]),"cutadapt5First"])))
  }
  #Check WT sequence is valid (ACGT characters only)
  all_characters <- unique(unlist(strsplit(dimsum_meta[["wildtypeSequence"]], "")))
  if(sum(!all_characters %in% c("A", "C", "G", "T"))!=0){
    stop("Invalid wild-type nucleotide sequence. Only valid nucleotide sequences allowed (A/C/T/G).", call. = FALSE)
  }
  return(exp_design)
}


