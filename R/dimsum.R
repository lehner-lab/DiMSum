
#' Run DiMSum pipeline
#'
#' This function runs the DiMSum pipeline.
#'
#' @param fastqFileDir Path to directory with input FASTQ files (required for WRAP)
#' @param fastqFileExtension FASTQ file extension (default:'.fastq')
#' @param gzipped Are FASTQ files are gzipped? (default:T)
#' @param stranded Is the library design stranded? (default:T)
#' @param paired Is the library design paired-end? (default:T)
#' @param barcodeDesignPath Path to barcode design file (tab-separated plain text file with barcode design)
#' @param barcodeErrorRate Maximum allowed error rate for the barcode (default:0.25)
#' @param experimentDesignPath Path to experimental design file (tab-separated plain text file with replicate structure)
#' @param experimentDesignPairDuplicates Are duplicate FASTQ files permitted in experimental design file? (default:F)
#' @param barcodeIdentityPath Path to barcode identity file (tab-separated plain text file mapping barcodes to variants)
#' @param countPath Path to variant count file for analysis with STEAM only (tab-separated plain text file with sample counts for all variants)
#' @param cutadaptCut5First cutadapt: remove bases from start of first read (before constant region trimming)
#' @param cutadaptCut5Second cutadapt: remove bases from start of second read (before constant region trimming)
#' @param cutadaptCut3First cutadapt: remove bases from end of first read (before constant region trimming)
#' @param cutadaptCut3Second cutadapt: remove bases from end of second read (before constant region trimming)
#' @param cutadapt5First cutadapt: sequence of 5' constant region to be trimmed (of the first read)
#' @param cutadapt5Second cutadapt: sequence of 5' constant region to be trimmed (of the second read)
#' @param cutadapt3First cutadapt: sequence of 3' constant region to be trimmed (of the first read; default: reverse complement of cutadapt5Second if not linked)
#' @param cutadapt3Second cutadapt: sequence of 3' constant region to be trimmed (of the second read; default: reverse complement of cutadapt5First if not linked)
#' @param cutadaptMinLength cutadapt: Discard reads shorter than LENGTH after trimming (default:50)
#' @param cutadaptErrorRate cutadapt: Maximum allowed error rate for trimming (default:0.2)
#' @param cutadaptOverlap cutadapt: Minimum overlap between read and constant region for trimming (default:3)
#' @param usearchMinQual USEARCH: minimum observed base quality to retain read pair (default:30)
#' @param usearchMaxee USEARCH: maximum number of expected errors to retain read pair (default:0.5)
#' @param usearchMinlen USEARCH: Discard pair if either read is shorter than this (default:64)
#' @param usearchMinovlen USEARCH: discard pair if alignment is shorter than given value (default:16)
#' @param outputPath Path to directory to use for output files
#' @param projectName Project name
#' @param wildtypeSequence Wild-type nucleotide sequence (A/C/G/T). Lower-case letters (a/c/g/t) indicate internal constant regions to be removed during WRAP.
#' @param permittedSequences A sequence of nucleotide codes (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N) with length matching the number of mutated positions i.e upper-case letters in wild-type nucleotide sequence (default:any base at mutated positions)
#' @param reverseComplement Reverse complement variant sequences before processing? (default:F)
#' @param sequenceType Coding potential of sequence; either noncoding/coding/auto (default:auto)
#' @param mutagenesisType Whether mutagenesis was performed at the nucleotide or codon/amino acid level; either random/codon (default:random)
#' @param transLibrary Trans library design i.e. read pairs correspond to distinct peptides with no overlap (default:F)
#' @param transLibraryReverseComplement Reverse complement second read in pair before concatenating (default:F)
#' @param bayesianDoubleFitness Improve double mutant fitness estimates using Bayesian framework (DISABLED: still in development)
#' @param bayesianDoubleFitnessLamD Poisson distribution for score likelihood (default:0.025)
#' @param fitnessMinInputCountAll Minimum input read count (in all replicate) to be retained during fitness calculations (default:0)
#' @param fitnessMinInputCountAny Minimum input read count (in any replicate) to be retained during fitness calculations (default:0)
#' @param fitnessMinOutputCountAll Minimum output read count (in all replicate) to be retained during fitness calculations (default:0)
#' @param fitnessMinOutputCountAny Minimum output read count (in any replicate) to be retained during fitness calculations (default:0)
#' @param fitnessHighConfidenceCount Minimum mean input read count for high confidence variants (default:10)
#' @param fitnessDoubleHighConfidenceCount Minimum input replicate read count for doubles used to derive prior for Bayesian doubles correction (default:50)
#' @param fitnessNormalise Normalise fitness values to minimise inter-replicate differences (default:T)
#' @param fitnessErrorModel Fit fitness error model (default:T)
#' @param maxSubstitutions Maximum number of nucleotide or amino acid substitutions for coding or non-coding sequences respectively (default:2)
#' @param mixedSubstitutions For coding sequences, are nonsynonymous variants with silent/synonymous substitutions in other codons allowed? (default:F)
#' @param retainIntermediateFiles Should intermediate files be retained? (default:F)
#' @param splitChunkSize FASTQ file split chunk size in bytes (default:3758096384)
#' @param retainedReplicates Comma-separated list of Input replicates (or experiment ids) to retain or 'all' (default:'all')
#' @param startStage Start at a specified pipeline stage (default:0)
#' @param stopStage Stop at a specified pipeline stage (default:5)
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
  experimentDesignPairDuplicates=F,
  barcodeIdentityPath,
  countPath,
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
  cutadaptOverlap=3,
  usearchMinQual=30,
  usearchMaxee=0.5,
  usearchMinlen=64,
  usearchMinovlen=16,
  outputPath,
  projectName,
  wildtypeSequence,
  permittedSequences,
  reverseComplement=F,
  sequenceType="auto",
  mutagenesisType="random",
  transLibrary=F,
  transLibraryReverseComplement=F,
  bayesianDoubleFitness=F,
  bayesianDoubleFitnessLamD=0.025,
  fitnessMinInputCountAll=0,
  fitnessMinInputCountAny=0,
  fitnessMinOutputCountAll=0,
  fitnessMinOutputCountAny=0,
  fitnessHighConfidenceCount=10,
  fitnessDoubleHighConfidenceCount=50,
  fitnessNormalise=T,
  fitnessErrorModel=T,
  maxSubstitutions=2,
  mixedSubstitutions=F,
  retainIntermediateFiles=F,
  splitChunkSize=3758096384,
  retainedReplicates="all",
  startStage=0,
  stopStage=5,
  numCores=1
  ){

  #Display welcome
  dimsum__status_message(paste("\n\n\n*******", "Running DiMSum pipeline", "*******\n\n\n"))
  dimsum__status_message(paste(formatDL(unlist(list(
    "Package version" = as.character(packageVersion("DiMSum")),
    "R version" = version$version.string))), collapse = "\n"))

  ### Basic checks
  ###########################

  #Only necessary if WRAP will be run
  if(is.null(countPath)){
    #Required binaries
    required_binaries <- c(
      "cutadapt", 
      "fastqc", 
      "gunzip", 
      "pandoc",
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
    dimsum__status_message(paste("\n\n\n*******", "Binary dependency versions", "*******\n\n\n"))
    dimsum__status_message(paste(formatDL(unlist(binary_versions)), collapse = "\n"))
  }

  ### Setup
  ###########################

  dimsum_arg_list <- list(
    "fastqFileDir" = list(fastqFileDir, c("character", "NULL")), #directory exists if running WRAP stages -- checked in dimsum__validate_input
    "fastqFileExtension" = list(fastqFileExtension, c("character")), #alphanumeric character string starting with '.' -- checked in dimsum__validate_input
    "gzipped" = list(gzipped, c("logical")), #logical -- checked in dimsum__validate_input
    "stranded" = list(stranded, c("logical")), #logical -- checked in dimsum__validate_input
    "paired" = list(paired, c("logical")), #logical -- checked in dimsum__validate_input
    "barcodeDesignPath" = list(barcodeDesignPath, c("character", "NULL")), #file exists (if not NULL)
    "barcodeErrorRate" = list(barcodeErrorRate, c("double")), #positive double less than 1 (zero inclusive) -- checked in dimsum__validate_input
    "experimentDesignPath" = list(experimentDesignPath, c("character")), #file exists -- checked in dimsum__get_experiment_design
    "experimentDesignPairDuplicates" = list(experimentDesignPairDuplicates, c("logical")), #logical -- checked in dimsum__validate_input
    "barcodeIdentityPath" = list(barcodeIdentityPath, c("character", "NULL")), #file exists (if not NULL) -- checked in dimsum__validate_input
    "countPath" = list(countPath, c("character", "NULL")), #file exists (if not NULL) -- checked in dimsum__validate_input
    "cutadaptCut5First" = list(cutadaptCut5First, c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadaptCut5Second" = list(cutadaptCut5Second, c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadaptCut3First" = list(cutadaptCut3First, c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadaptCut3Second" = list(cutadaptCut3Second, c("integer", "NULL")), #strictly positive integer (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadapt5First" = list(cutadapt5First, c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadapt5Second" = list(cutadapt5Second, c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadapt3First" = list(cutadapt3First, c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadapt3Second" = list(cutadapt3Second, c("character", "NULL")), #AGCT character string (if not NULL) -- checked in dimsum__check_experiment_design
    "cutadaptMinLength" = list(cutadaptMinLength, c("integer")), #strictly positive integer -- checked in dimsum__check_experiment_design
    "cutadaptErrorRate" = list(cutadaptErrorRate, c("double")), #positive double less than 1 (zero inclusive) -- checked in dimsum__check_experiment_design
    "cutadaptOverlap" = list(cutadaptOverlap, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__check_experiment_design
    "usearchMinQual" = list(usearchMinQual, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "usearchMaxee" = list(usearchMaxee, c("double")), #strictly positive double -- checked in dimsum__validate_input
    "usearchMinlen" = list(usearchMinlen, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "usearchMinovlen" = list(usearchMinovlen, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "outputPath" = list(outputPath, c("character")), #directory exists -- checked in dimsum__validate_input
    "projectName" = list(projectName, c("character")), #character string -- checked in dimsum__validate_input
    "wildtypeSequence" = list(wildtypeSequence, c("character")), #AaGgCcTt character string -- checked in dimsum__validate_input
    "permittedSequences" = list(permittedSequences, c("character", "NULL")), #ACGTRYSWKMBDHVN character string -- checked in dimsum__validate_input
    "reverseComplement" = list(reverseComplement, c("logical")), #logical -- checked in dimsum__validate_input
    "sequenceType" = list(sequenceType, c("character")), #character string; either noncoding/coding/auto -- checked in dimsum__validate_input
    "mutagenesisType" = list(mutagenesisType, c("character")), #character string; either random/codon -- checked in dimsum__validate_input
    "transLibrary" = list(transLibrary, c("logical")), #logical -- checked in dimsum__validate_input
    "transLibraryReverseComplement" = list(transLibraryReverseComplement, c("logical")), #logical -- checked in dimsum__validate_input
    "bayesianDoubleFitness" = list(bayesianDoubleFitness, c("logical")), #logical -- checked in dimsum__validate_input
    "bayesianDoubleFitnessLamD" = list(bayesianDoubleFitnessLamD, c("double")), #strictly positive double -- checked in dimsum__validate_input
    "fitnessMinInputCountAll" = list(fitnessMinInputCountAll, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinInputCountAny" = list(fitnessMinInputCountAny, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinOutputCountAll" = list(fitnessMinOutputCountAll, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinOutputCountAny" = list(fitnessMinOutputCountAny, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessHighConfidenceCount" = list(fitnessHighConfidenceCount, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessDoubleHighConfidenceCount" = list(fitnessDoubleHighConfidenceCount, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessNormalise" = list(fitnessNormalise, c("logical")), #logical -- checked in dimsum__validate_input
    "fitnessErrorModel" = list(fitnessErrorModel, c("logical")), #logical -- checked in dimsum__validate_input
    "maxSubstitutions" = list(maxSubstitutions, c("integer")), #positive integer (greater than 1) -- checked in dimsum__validate_input
    "mixedSubstitutions" = list(mixedSubstitutions, c("logical")), #logical -- checked in dimsum__validate_input
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

  #Initialise reports directory (and erase if already exists)
  dimsum__create_dir(file.path(exp_metadata[["project_path"]], "reports"), execute = T)
  #Render report (preemptively)
  dimsum__render_report(dimsum_meta = exp_metadata, initialise = TRUE)

  ### Pipeline stages
  ###########################

  ### Start pipeline tracking
  pipeline <- list()
  pipeline[['initial']] <- exp_metadata

  ### Stage 0 (WRAP): Run demultiplex on all fastq files
  pipeline[['0_demultiplex']] <- dimsum_stage_demultiplex(dimsum_meta = pipeline[['initial']], demultiplex_outpath = file.path(pipeline[['initial']][["tmp_path"]], "0_demultiplex"))

  ### Stage 1.1 (WRAP): Run FASTQC on all fastq files
  pipeline[['1_qualitycontrol']] <- dimsum_stage_fastqc(dimsum_meta = pipeline[['0_demultiplex']], fastqc_outpath = file.path(pipeline[['0_demultiplex']][["tmp_path"]], "1_qualitycontrol"), 
    report_outpath = file.path(pipeline[['0_demultiplex']][["project_path"]], "reports"))

  ### Stage 1.2 (WRAP): Unzip FASTQ files if necessary
  pipeline[['1_unzip']] <- dimsum_stage_unzip(dimsum_meta = pipeline[['1_qualitycontrol']], fastq_outpath = file.path(pipeline[['1_qualitycontrol']][["tmp_path"]], "1_unzip"))

  ### Stage 1.3 (WRAP): Split FASTQ files
  pipeline[['1_split']] <- dimsum_stage_split(dimsum_meta = pipeline[['1_unzip']], split_outpath = file.path(pipeline[['1_unzip']][["tmp_path"]], "1_split"))

  ### Stage 2 (WRAP): Remove constant regions from FASTQ files with cutadapt if necessary
  pipeline[['2_trim']] <- dimsum_stage_cutadapt(dimsum_meta = pipeline[['1_split']], cutadapt_outpath = file.path(pipeline[['1_split']][["tmp_path"]], "2_trim"), 
    report_outpath = file.path(pipeline[['1_split']][["project_path"]], "reports"))

  ### Stage 3.1 (WRAP): Merge paired-end reads with USEARCH
  pipeline[['3_align']] <- dimsum_stage_usearch(dimsum_meta = pipeline[['2_trim']], usearch_outpath = file.path(pipeline[['2_trim']][["tmp_path"]], "3_align"), 
    report_outpath = file.path(pipeline[['2_trim']][["project_path"]], "reports"))

  ### Stage 3.2 (WRAP): Tally unique read counts with FASTX-Toolkit
  pipeline[['3_tally']] <- dimsum_stage_unique(dimsum_meta = pipeline[['3_align']], unique_outpath = file.path(pipeline[['3_align']][["tmp_path"]], "3_tally"))

  ### Stage 4 (STEAM): Merge count files, process variant sequence and filter for desired substitution variants
  pipeline[['4_process']] <- dimsum_stage_merge(dimsum_meta = pipeline[['3_tally']], merge_outpath = pipeline[['3_tally']][["project_path"]], 
    report_outpath = file.path(pipeline[['3_tally']][["project_path"]], "reports"))

  ### Stage 5 (STEAM): Calculate fitness and model error
  pipeline[['5_analyse']] <- dimsum_stage_counts_to_fitness(dimsum_meta = pipeline[['4_process']], fitness_outpath = pipeline[['4_process']][["project_path"]], 
    report_outpath = file.path(pipeline[['4_process']][["project_path"]], "reports"))

  ### Save workspace
  ###########################

  dimsum__status_message("\n\n\nSaving workspace image...\n")
  dimsum__save_metadata(dimsum_meta = pipeline[['5_analyse']], n = 1)
  dimsum__render_report(dimsum_meta = pipeline[['5_analyse']], finalise = TRUE)
  dimsum__status_message("Done\n")

}


