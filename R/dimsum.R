
#' Run DiMSum pipeline
#'
#' This function runs the DiMSum pipeline.
#'
#' @param runDemo Run the DiMSum demo (default:F)
#' @param fastqFileDir Path to directory containing input FASTQ files (required for WRAP)
#' @param fastqFileExtension FASTQ file extension (default:'.fastq')
#' @param gzipped Are FASTQ files are gzipped? (default:T)
#' @param stranded Is the library design stranded? (default:T)
#' @param paired Is the library design paired-end? (default:T)
#' @param barcodeDesignPath Path to barcode design file (tab-separated plain text file with barcode design)
#' @param barcodeErrorRate Maximum allowed error rate for barcode to be matched (default:0.25)
#' @param experimentDesignPath Path to Experimental Design File (required if '--runDemo'=F)
#' @param experimentDesignPairDuplicates Are multiple instances of FASTQ files in the Experimental Design File permitted? (default:F)
#' @param barcodeIdentityPath Path to Variant Identity File (tab-separated plain text file mapping barcodes to variants)
#' @param countPath Path to Variant Count File for analysis with STEAM only (tab-separated plain text file with sample counts for all variants)
#' @param synonymSequencePath Path to Synonym Sequences File with coding sequences for which synonymous variant fitness should be quantified (default: plain text file with one coding nucleotide sequence per line)
#' @param cutadaptCut5First Remove fixed number of bases from start (5') of first (or only) read before constant region trimming (optional)
#' @param cutadaptCut5Second Remove fixed number of bases from start (5') of second read in pair before constant region trimming (optional)
#' @param cutadaptCut3First Remove fixed number of bases from end (3') of first (or only) read before constant region trimming (optional)
#' @param cutadaptCut3Second Remove fixed number of bases from end (3') of second read in pair before constant region trimming (optional)
#' @param cutadapt5First Sequence of 5' constant region to be trimmed from first (or only) read (optional)
#' @param cutadapt5Second Sequence of 5' constant region to be trimmed from second read in pair (optional)
#' @param cutadapt3First Sequence of 3' constant region to be trimmed from first (or only) read (default: reverse complement of '--cutadapt5Second')
#' @param cutadapt3Second Sequence of 3' constant region to be trimmed from second read in pair (default: reverse complement of '--cutadapt5First')
#' @param cutadaptMinLength Discard reads shorter than LENGTH after trimming (default:50)
#' @param cutadaptErrorRate Maximum allowed error rate for trimming constant regions (default:0.2)
#' @param cutadaptOverlap Minimum overlap between read and constant region for trimming (default:3)
#' @param vsearchMinQual Minimum Phred base quality score required to retain read or read pair (default:30)
#' @param vsearchMaxee Maximum number of expected errors tolerated to retain read or read pair (default:0.5)
#' @param vsearchMinovlen Discard read pair if the alignment length is shorter than this (default:10)
#' @param outputPath Path to directory to use for output files (default:'./' i.e. current working directory)
#' @param projectName Project name and directory where results are to be saved (default:'DiMSum_Project')
#' @param wildtypeSequence Wild-type nucleotide sequence (A/C/G/T). Lower-case bases (a/c/g/t) indicate internal constant regions to be removed (required if '--runDemo'=F)
#' @param permittedSequences Nucleotide sequence of IUPAC ambiguity codes (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N) with length matching the number of mutated positions (i.e upper-case letters) in '--wildtypeSequence' (default:N i.e. any substitution mutation allowed)
#' @param reverseComplement Reverse complement sequence (default:F)
#' @param sequenceType Coding potential of sequence: either 'noncoding', 'coding' or 'auto'. If the specified wild-type nucleotide sequence ('--wildtypeSequence') has a valid translation without a premature STOP codon, it is assumed to be 'coding' (default:'auto')
#' @param mutagenesisType Whether mutagenesis was performed at the nucleotide or codon/amino acid level; either 'random' or 'codon' (default:'random')
#' @param transLibrary Paired-end reads correspond to distinct molecules? (default:F)
#' @param transLibraryReverseComplement Reverse complement second read in pair (default:F)
#' @param bayesianDoubleFitness In development: improve double mutant fitness estimates using Bayesian framework (DISABLED: still in development)
#' @param bayesianDoubleFitnessLamD In development: Poisson distribution for score likelihood (default:0.025)
#' @param fitnessMinInputCountAll Minimum input read count (in all replicates) to be retained during fitness calculations (default:0)
#' @param fitnessMinInputCountAny Minimum input read count (in any replicate) to be retained during fitness calculations (default:0)
#' @param fitnessMinOutputCountAll Minimum output read count (in all replicates) to be retained during fitness calculations (default:0)
#' @param fitnessMinOutputCountAny Minimum output read count (in any replicates) to be retained during fitness calculations (default:0)
#' @param fitnessHighConfidenceCount In development: minimum mean input read count for high confidence variants (default:10)
#' @param fitnessDoubleHighConfidenceCount In development: minimum input replicate read count for doubles used to derive prior for Bayesian doubles correction (default:50)
#' @param fitnessNormalise Normalise fitness values to minimise inter-replicate differences (default:T)
#' @param fitnessErrorModel Fit fitness error model (default:T)
#' @param indels Indel variants to be retained: either 'all', 'none' or a comma-separated list of sequence lengths (default:'none')
#' @param maxSubstitutions Maximum number of nucleotide or amino acid substitutions for coding or non-coding sequences respectively (default:2)
#' @param mixedSubstitutions For coding sequences, are nonsynonymous variants with silent/synonymous substitutions in other codons allowed? (default:F)
#' @param retainIntermediateFiles Should intermediate files be retained? Intermediate files can be many gigabytes, but are required to rerun DiMSum starting at intermediate pipeline stages (default:F)
#' @param splitChunkSize Internal: FASTQ file split chunk size in bytes (default:3758096384)
#' @param retainedReplicates Comma-separated list of (integer) experiment replicates to retain or 'all' (default:'all')
#' @param startStage (Re-)Start DiMSum at a specific pipeline stage (default:0)
#' @param stopStage Stop DiMSum at a specific pipeline stage (default:5)
#' @param numCores Number of available CPU cores. All pipeline stages make use of parallel computing to decrease runtime if multiple cores are available (default:1)
#'
#' @return Nothing
#' @export
dimsum <- function(
  runDemo=F,
  fastqFileDir=NULL,
  fastqFileExtension=".fastq",
  gzipped=T,
  stranded=T,
  paired=T,
  barcodeDesignPath=NULL,
  barcodeErrorRate=0.25,
  experimentDesignPath=NULL,
  experimentDesignPairDuplicates=F,
  barcodeIdentityPath=NULL,
  countPath=NULL,
  synonymSequencePath=NULL,
  cutadaptCut5First=NULL,
  cutadaptCut5Second=NULL,
  cutadaptCut3First=NULL,
  cutadaptCut3Second=NULL,
  cutadapt5First=NULL,
  cutadapt5Second=NULL,
  cutadapt3First=NULL,
  cutadapt3Second=NULL,
  cutadaptMinLength=50,
  cutadaptErrorRate=0.2,
  cutadaptOverlap=3,
  vsearchMinQual=30,
  vsearchMaxee=0.5,
  vsearchMinovlen=10,
  outputPath="./",
  projectName="DiMSum_Project",
  wildtypeSequence=NULL,
  permittedSequences=NULL,
  reverseComplement=F,
  sequenceType="auto",
  mutagenesisType="random",
  transLibrary=F,
  transLibraryReverseComplement=F,
  bayesianDoubleFitness=F,
  bayesianDoubleFitnessLamD=0.025,
  fitnessMinInputCountAll="0",
  fitnessMinInputCountAny="0",
  fitnessMinOutputCountAll="0",
  fitnessMinOutputCountAny="0",
  fitnessHighConfidenceCount=10,
  fitnessDoubleHighConfidenceCount=50,
  fitnessNormalise=T,
  fitnessErrorModel=T,
  indels='none',
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
    "R version" = version$version.string))), collapse = "\n"), newline = T)

  ### Demo
  ###########################

  if(runDemo){
    wildtypeSequence <- "GGTAATAGCAGAGGGGGTGGAGCTGGTTTGGGAAACAATCAAGGTAGTAATATGGGTGGTGGGATGAACTTTGGTGCGTTCAGCATTAATCCAGCCATGATGGCTGCCGCCCAGGCAGCACTACAG"
    cutadapt5First="TGGCTTTGGGAATCAGGGTGGATTT"
    cutadapt5Second="ACATGCCCATCATACCCCAACTGCT"
    experimentDesignPath <- system.file("demo", "experimentDesign_Toy.txt", package = "DiMSum")
    if(is.null(fastqFileDir)){
      countPath <- system.file("demo", "countFile_Toy.txt", package = "DiMSum")
    }
  }

  ### Basic checks
  ###########################

  #If WRAP will be run
  if(is.null(countPath)){
    #Required binaries
    required_binaries <- c(
      "cutadapt", 
      "fastqc", 
      "gunzip", 
      "pandoc",
      "vsearch",
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
      Pandoc = system("pandoc -v", intern = TRUE)[1],      
      VSEARCH = unlist(strsplit(system("vsearch --version 2>&1", intern = TRUE), ","))[1],
      starcode = system("starcode --version 2>&1", intern = TRUE))))

    #Display binary versions
    dimsum__status_message(paste("\n\n\n*******", "Binary dependency versions", "*******\n\n\n"))
    dimsum__status_message(paste(formatDL(unlist(binary_versions)), collapse = "\n"), newline = T)
  }

  #If WRAP will NOT be run
  if(!is.null(countPath)){
    #Required binaries
    required_binaries <- c(
      "pandoc")
    which_binaries <- Sys.which(required_binaries)
    missing_binaries <- names(which_binaries)[which_binaries==""]
    if(length(missing_binaries)!=0){
      stop(paste0("Required executables not installed. Please install the following software: ", paste(missing_binaries, sep = ", ")), call. = FALSE)
    }

    #Binary versions
    suppressMessages(suppressWarnings(binary_versions <- list(
      Pandoc = system("pandoc -v", intern = TRUE)[1])))

    #Display binary versions
    dimsum__status_message(paste("\n\n\n*******", "Binary dependency versions", "*******\n\n\n"))
    dimsum__status_message(paste(formatDL(unlist(binary_versions)), collapse = "\n"), newline = T)
  }

  ### Setup
  ###########################

  dimsum_arg_list <- list(
    "runDemo" = list(runDemo, c("logical")), #logical -- checked in dimsum__validate_input
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
    "synonymSequencePath" = list(synonymSequencePath, c("character", "NULL")), #file exists (if not NULL) -- checked in dimsum__validate_input
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
    "vsearchMinQual" = list(vsearchMinQual, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
    "vsearchMaxee" = list(vsearchMaxee, c("double")), #strictly positive double -- checked in dimsum__validate_input
    "vsearchMinovlen" = list(vsearchMinovlen, c("integer")), #strictly positive integer -- checked in dimsum__validate_input
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
    "fitnessMinInputCountAll" = list(fitnessMinInputCountAll, c("character")), #character string; positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinInputCountAny" = list(fitnessMinInputCountAny, c("character")), #character string; positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinOutputCountAll" = list(fitnessMinOutputCountAll, c("character")), #character string; positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessMinOutputCountAny" = list(fitnessMinOutputCountAny, c("character")), #character string; positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessHighConfidenceCount" = list(fitnessHighConfidenceCount, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessDoubleHighConfidenceCount" = list(fitnessDoubleHighConfidenceCount, c("integer")), #positive integer (zero inclusive) -- checked in dimsum__validate_input
    "fitnessNormalise" = list(fitnessNormalise, c("logical")), #logical -- checked in dimsum__validate_input
    "fitnessErrorModel" = list(fitnessErrorModel, c("logical")), #logical -- checked in dimsum__validate_input
    "indels" = list(indels, c("character")), #character string; either all/none or a comma-separated list of sequence lengths -- checked in dimsum__validate_input
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

  #Check count file
  dimsum__check_countfile(exp_metadata)

  #Check barcode identity file
  dimsum__check_barcodeidentityfile(exp_metadata)

  #Get and check synonym sequences
  exp_metadata[["synonym_sequences"]] <- dimsum__get_synonymsequences(exp_metadata)

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

  ### Stage 3.1 (WRAP): Merge paired-end reads with VSEARCH
  pipeline[['3_align']] <- dimsum_stage_vsearch(dimsum_meta = pipeline[['2_trim']], vsearch_outpath = file.path(pipeline[['2_trim']][["tmp_path"]], "3_align"), 
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


