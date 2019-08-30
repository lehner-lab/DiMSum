
#' Run DiMSum pipeline
#'
#' This function runs the DiMSum pipeline.
#'
#' @param fastqFileDir Path to directory with input FASTQ files
#' @param fastqFileExtension FASTQ file extension (default:".fastq")
#' @param gzipped Are FASTQ files are gzipped? (default:T)
#' @param stranded Is the library design stranded? (default:T)
#' @param paired Is the library design paired-end? (default:T)
#' @param barcodeDesignPath Path to barcode design file (tab-separated plain text file with barcode design)
#' @param barcodeErrorRate Maximum allowed error rate for the barcode (default:0.25)
#' @param experimentDesignPath Path to experimental design file (tab-separated plain text file with replicate structure)
#' @param cutadaptCut5First cutadapt: remove bases from start of first read (before constant region trimming)
#' @param cutadaptCut5Second cutadapt: remove bases from start of second read (before constant region trimming)
#' @param cutadaptCut3First cutadapt: remove bases from end of first read (before constant region trimming)
#' @param cutadaptCut3Second cutadapt: remove bases from end of second read (before constant region trimming)
#' @param cutadapt5First cutadapt: sequence of 5' constant region (of the first read)
#' @param cutadapt5Second cutadapt: sequence of 5' constant region (of the second read)
#' @param cutadapt3First cutadapt: sequence of 3' constant region (of the first read; default: reverse complement of cutadapt5Second)
#' @param cutadapt3Second cutadapt: sequence of 3' constant region (of the second read; default: reverse complement of cutadapt5First)
#' @param cutadaptMinLength cutadapt: Discard reads shorter than LENGTH (default:50)
#' @param cutadaptErrorRate cutadapt: Maximum allowed error rate (default:0.2)
#' @param usearchMinQual USEARCH: minimum observed base quality to retain read pair
#' @param usearchMaxee USEARCH: maximum number of expected errors to retain read pair
#' @param usearchMinlen USEARCH: Discard pair if either read is shorter than this (default:64)
#' @param usearchMinovlen USEARCH: discard pair if alignment is shorter than given value (default:16)
#' @param outputPath Path to directory to use for output files
#' @param projectName Project name
#' @param wildtypeSequence Wild-type nucleotide sequence (A/C/G/T). Lower-case letters (a/c/g/t) indicate internal constant regions to be removed before fitness calculations.
#' @param sequenceType Coding potential of sequence; either noncoding/coding/auto (default:auto)
#' @param transLibrary Trans library design i.e. read pairs correspond to distinct peptides (no overlap)
#' @param bayesianDoubleFitness Improve double mutant fitness estimates using Bayesian framework (default:F)
#' @param bayesianDoubleFitnessLamD Poisson distribution for score likelihood (default:0.025)
#' @param fitnessMinInputCountAll Minimum input read count (in all replicate) to be retained during fitness calculations (default:0)
#' @param fitnessMinInputCountAny Minimum input read count (in any replicate) to be retained during fitness calculations (default:5)
#' @param fitnessHighConfidenceCount Minimum mean input read count for high confidence variants (default:10)
#' @param fitnessDoubleHighConfidenceCount Minimum input replicate read count for doubles used to derive prior for Bayesian doubles correction (default:50)
#' @param fitnessMaxSubstitutions Maximum number of nucleotide or amino acid substitutions for coding or non-coding sequences respectively (default:2)
#' @param retainIntermediateFiles Should intermediate files be retained? (default:F)
#' @param splitChunkSize FASTQ file split chunk size in bytes (default:3758096384)
#' @param retainedReplicates Comma-separated list of Input replicates (or experiment ids) to retain or 'all' (default:'all')
#' @param startStage Start at a specified pipeline stage (default:1)
#' @param stopStage Stop at a specified pipeline stage (default:0 i.e. no stop condition)
#' @param numCores Number of available CPU cores (default:1)
#'
#' @return Nothing
#' @export
dimsum <- function(
  fastqFileDir,
  fastqFileExtension=".fastq",
  gzipped=T,
  stranded=T,
  paired=T,
  barcodeDesignPath,
  barcodeErrorRate=0.25,
  experimentDesignPath,
  cutadaptCut5First,
  cutadaptCut5Second,
  cutadaptCut3First,
  cutadaptCut3Second,
  cutadapt5First,
  cutadapt5Second,
  cutadapt3First,
  cutadapt3Second,
  cutadaptMinLength=50,
  cutadaptErrorRate=0.2,
  usearchMinQual,
  usearchMaxee,
  usearchMinlen=64,
  usearchMinovlen=16,
  outputPath,
  projectName,
  wildtypeSequence,
  sequenceType="auto",
  transLibrary=F,
  bayesianDoubleFitness=F,
  bayesianDoubleFitnessLamD=0.025,
  fitnessMinInputCountAll=0,
  fitnessMinInputCountAny=5,
  fitnessHighConfidenceCount=10,
  fitnessDoubleHighConfidenceCount=50,
  fitnessMaxSubstitutions=2,
  retainIntermediateFiles=F,
  splitChunkSize=3758096384,
  retainedReplicates="all",
  startStage=1,
  stopStage=0,
  numCores=1
  ){

  #Display welcome
  message(paste("\n\n\n*******", "Running DiMSum pipeline", "*******\n\n\n"))
  message(paste(formatDL(unlist(list("Package version" = as.character(packageVersion("DiMSum"))))), collapse = "\n"))
  message(paste(formatDL(unlist(list("R version" = version$version.string))), collapse = "\n"))

  ### Basic checks
  ###########################

  #Required binaries
  required_binaries <- c(
    "cat", 
    "cp", 
    "cutadapt", 
    "fastqc", 
    "gunzip", 
    "head", 
    "usearch",
    "starcode")
  which_binaries <- Sys.which(required_binaries)
  missing_binaries <- names(which_binaries)[which_binaries==""]
  if(length(missing_binaries)!=0){
    stop(paste0("Required executables not installed. Please install the following software: ", paste(missing_binaries, sep = ", ")), call. = FALSE)
  }

  #Binary versions
  suppressMessages(suppressWarnings(binary_versions <- list(
    cutadapt = rev(system("cutadapt --version", intern = TRUE))[1],
    FastQC = rev(system("fastqc --version", intern = TRUE))[1],
    USEARCH = rev(system("usearch --version", intern = TRUE))[1],
    starcode = system("starcode --version 2>&1", intern = TRUE))))

  #Display binary versions
  message(paste("\n\n\n*******", "Binary dependency versions", "*******\n\n\n"))
  message(paste(formatDL(unlist(binary_versions)), collapse = "\n"))

  ### Setup
  ###########################

  dimsum_arg_list <- list(
    "fastqFileDir" = list(fastqFileDir, c("character")), #directory exists -- checked in dimsum__validate_input
    "fastqFileExtension" = list(fastqFileExtension, c("character")), #alphanumeric character string starting with '.' -- checked in dimsum__validate_input
    "gzipped" = list(gzipped, c("logical")), #logical -- checked in dimsum__validate_input
    "stranded" = list(stranded, c("logical")), #logical -- checked in dimsum__validate_input
    "paired" = list(paired, c("logical")), #logical -- checked in dimsum__validate_input
    "barcodeDesignPath" = list(barcodeDesignPath, c("character", "NULL")), #file exists (if not NULL)
    "barcodeErrorRate" = list(barcodeErrorRate, c("double")), #positive double (zero inclusive)
    "experimentDesignPath" = list(experimentDesignPath, c("character")), #file exists -- checked in dimsum__get_experiment_design
    "cutadaptCut5First" = list(cutadaptCut5First, c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__get_experiment_design
    "cutadaptCut5Second" = list(cutadaptCut5Second, c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__get_experiment_design
    "cutadaptCut3First" = list(cutadaptCut3First, c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__get_experiment_design
    "cutadaptCut3Second" = list(cutadaptCut3Second, c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__get_experiment_design
    "cutadapt5First" = list(cutadapt5First, c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__get_experiment_design
    "cutadapt5Second" = list(cutadapt5Second, c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__get_experiment_design
    "cutadapt3First" = list(cutadapt3First, c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__get_experiment_design
    "cutadapt3Second" = list(cutadapt3Second, c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__get_experiment_design
    "cutadaptMinLength" = list(cutadaptMinLength, c("integer")), #strictly positive integer -- checked in dimsum__get_experiment_design
    "cutadaptErrorRate" = list(cutadaptErrorRate, c("double")), #positive double (zero inclusive) -- checked in dimsum__get_experiment_design
    "usearchMinQual" = list(usearchMinQual, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "usearchMaxee" = list(usearchMaxee, c("double")), #strictly positive double -- checked in dimsum__validate_input
    "usearchMinlen" = list(usearchMinlen, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "usearchMinovlen" = list(usearchMinovlen, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "outputPath" = list(outputPath, c("character")), #directory exists -- checked in dimsum__validate_input
    "projectName" = list(projectName, c("character")), #character string -- checked in dimsum__validate_input
    "wildtypeSequence" = list(wildtypeSequence, c("character")), #AGCT character string -- checked in dimsum__validate_input
    "sequenceType" = list(sequenceType, c("character")), #character string; either noncoding/coding/auto -- checked in dimsum__validate_input
    "transLibrary" = list(transLibrary, c("logical")), #logical -- checked in dimsum__validate_input
    "bayesianDoubleFitness" = list(bayesianDoubleFitness, c("logical")), #logical -- checked in dimsum__validate_input
    "bayesianDoubleFitnessLamD" = list(bayesianDoubleFitnessLamD, c("double")), #strictly positive double -- checked in dimsum__validate_input
    "fitnessMinInputCountAll" = list(fitnessMinInputCountAll, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinInputCountAny" = list(fitnessMinInputCountAny, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessHighConfidenceCount" = list(fitnessHighConfidenceCount, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessDoubleHighConfidenceCount" = list(fitnessDoubleHighConfidenceCount, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMaxSubstitutions" = list(fitnessMaxSubstitutions, c("integer")), #positive integer (greater than 1) -- checked in dimsum__validate_input
    "retainIntermediateFiles" = list(retainIntermediateFiles, c("logical")), #logical -- checked in dimsum__validate_input
    "splitChunkSize" = list(splitChunkSize, c("double")), #strictly positive double -- checked in dimsum__validate_input
    "retainedReplicates" = list(retainedReplicates, c("character")), #comma-separated list of integers or "all" -- checked in dimsum__validate_input
    "startStage" = list(startStage, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "stopStage" = list(stopStage, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "numCores" = list(numCores, c("integer")) #strictly positive integer -- checked in dimsum__validate_input
    )

  #Validate input
  exp_metadata <- dimsum__validate_input(dimsum_arg_list)

  #Get and check barcode design (if provided)
  exp_metadata[["barcode_design"]] <- dimsum__get_barcode_design(exp_metadata)

  #Get and check experiment design
  exp_metadata[["exp_design"]] <- dimsum__get_experiment_design(exp_metadata)

  #Create project directory (if doesn't already exist)
  exp_metadata[["project_path"]] <- file.path(exp_metadata[["outputPath"]], exp_metadata[["projectName"]])
  suppressWarnings(dir.create(exp_metadata[["project_path"]]))
  #Create temp directory (if doesn't already exist)
  exp_metadata[["tmp_path"]] <- file.path(exp_metadata[["project_path"]], "tmp")
  suppressWarnings(dir.create(exp_metadata[["tmp_path"]]))

  #inisiDelete intermediate files string
  exp_metadata[["deleteIntermediateFiles"]] <- NULL

  ### Pipeline stages
  ###########################

  ### Step 0: Start pipeline tracking
  pipeline <- list()
  pipeline[['0_original']] <- exp_metadata

  ### Step 1: Run demultiplex on all fastq files
  pipeline[['1_demultiplex']] <- dimsum_stage_demultiplex(dimsum_meta = pipeline[['0_original']], demultiplex_outpath = file.path(pipeline[['0_original']][["tmp_path"]], "demultiplex"))

  ### Step 2.1: Run FASTQC on all fastq files
  pipeline[['2_fastqc']] <- dimsum_stage_fastqc(dimsum_meta = pipeline[['1_demultiplex']], fastqc_outpath = file.path(pipeline[['1_demultiplex']][["tmp_path"]], "fastqc"), 
    report_outpath = file.path(pipeline[['1_demultiplex']][["project_path"]], "reports"))

  ### Step 2.2: Unzip FASTQ files if necessary
  pipeline[['2_fastq']] <- dimsum_stage_unzip(dimsum_meta = pipeline[['2_fastqc']], fastq_outpath = file.path(pipeline[['2_fastqc']][["tmp_path"]], "fastq"))

  ### Step 2.3: Split FASTQ files
  pipeline[['2_split']] <- dimsum_stage_split(dimsum_meta = pipeline[['2_fastq']], split_outpath = file.path(pipeline[['2_fastq']][["tmp_path"]], "split"))

  ### Step 3: Remove adapters from FASTQ files with cutadapt if necessary
  pipeline[['3_cutadapt']] <- dimsum_stage_cutadapt(dimsum_meta = pipeline[['2_split']], cutadapt_outpath = file.path(pipeline[['2_split']][["tmp_path"]], "cutadapt"), 
    report_outpath = file.path(pipeline[['2_split']][["project_path"]], "reports"))

  ### Step 4: Merge paired-end reads with USEARCH
  pipeline[['4_usearch']] <- dimsum_stage_usearch(dimsum_meta = pipeline[['3_cutadapt']], usearch_outpath = file.path(pipeline[['3_cutadapt']][["tmp_path"]], "usearch"), 
    report_outpath = file.path(pipeline[['3_cutadapt']][["project_path"]], "reports"))

  ### Step 5: Get unique aligned read counts with FASTX-Toolkit
  pipeline[['5_unique']] <- dimsum_stage_unique(dimsum_meta = pipeline[['4_usearch']], unique_outpath = file.path(pipeline[['4_usearch']][["tmp_path"]], "unique"))

  ### Step 6: Merge variant count tables
  pipeline[['6_merge']] <- dimsum_stage_merge(dimsum_meta = pipeline[['5_unique']], merge_outpath = pipeline[['5_unique']][["project_path"]], 
    report_outpath = file.path(pipeline[['5_unique']][["project_path"]], "reports"))

  ### Step 7: Calculate fitness
  pipeline[['7_fitness']] <- dimsum_stage_counts_to_fitness(dimsum_meta = pipeline[['6_merge']], fitness_outpath = pipeline[['6_merge']][["project_path"]], 
    report_outpath = file.path(pipeline[['6_merge']][["project_path"]], "reports"))

  ### Save workspace
  ###########################

  message("\n\n\nSaving workspace image...")
  dimsum__save_metadata(dimsum_meta = pipeline[['7_fitness']], n = 1)
  message("Done")

  ### Save report html
  ###########################

  message("\n\n\nSaving summary report...")
  write(dimsum__reports_summary(dimsum_meta = pipeline[['7_fitness']]), file = file.path(pipeline[['7_fitness']][["project_path"]], "reports_summary.html"))
  message("Done")
}


