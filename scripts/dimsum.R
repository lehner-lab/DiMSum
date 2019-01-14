
#' Run DiMSum pipeline
#'
#' This function runs the DiMSum pipeline.
#'
#' @param fastqFileDir Path to directory with input FASTQ files
#' @param fastqFileExtension FASTQ file extension
#' @param gzipped Are FASTQ files are gzipped?
#' @param stranded Is the library design stranded?
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
#' @param usearchAttemptExactMinovlen USEARCH: Attempt exact alignment of --usearchMinovlen (default:F)
#' @param outputPath Path to directory to use for output files
#' @param projectName Project name
#' @param wildtypeSequence Wild-type nucleotide sequence
#' @param transLibrary Trans library design i.e. read pairs correspond to distinct peptides (no overlap)
#' @param startStage Start at a specified pipeline stage
#' @param stopStage "Stop at specified pipeline stage (default:0, no stop condition)
#' @param numCores Number of available CPU cores
#'
#' @return Nothing
#' @export
dimsum <- function(
  fastqFileDir,
  fastqFileExtension=".fastq",
  gzipped=T,
  stranded=T,
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
  usearchAttemptExactMinovlen=F,
  outputPath,
  projectName,
  wildtypeSequence,
  transLibrary=F,
  startStage=1,
  stopStage=0,
  numCores=1
  ){

  #Metadata object
  exp_metadata <- list()

  ### Save metadata
  #Remove trailing "/" if present
  exp_metadata[["fastq_path_original"]] <- gsub("/$", "", fastqFileDir)
  exp_metadata[["fastq_file_extension"]] <- fastqFileExtension
  exp_metadata[["gzipped"]] <- gzipped
  exp_metadata[["stranded"]] <- stranded
  exp_metadata[["barcode_design_path"]] <- barcodeDesignPath
  exp_metadata[["barcodeErrorRate"]] <- barcodeErrorRate
  exp_metadata[["experiment_design_path"]] <- experimentDesignPath
  exp_metadata[["cutadaptCut5First"]] <- cutadaptCut5First
  exp_metadata[["cutadaptCut5Second"]] <- cutadaptCut5Second
  exp_metadata[["cutadaptCut3First"]] <- cutadaptCut3First
  exp_metadata[["cutadaptCut3Second"]] <- cutadaptCut3Second
  exp_metadata[["cutadapt5First"]] <- cutadapt5First
  exp_metadata[["cutadapt5Second"]] <- cutadapt5Second
  exp_metadata[["cutadapt3First"]] <- cutadapt3First
  exp_metadata[["cutadapt3Second"]] <- cutadapt3Second
  exp_metadata[["cutadaptMinLength"]] <- cutadaptMinLength
  exp_metadata[["cutadaptErrorRate"]] <- cutadaptErrorRate
  exp_metadata[["usearchMinQual"]] <- usearchMinQual
  exp_metadata[["usearchMaxee"]] <- usearchMaxee
  exp_metadata[["usearchMinlen"]] <- usearchMinlen
  exp_metadata[["usearchMinovlen"]] <- usearchMinovlen
  exp_metadata[["usearchAttemptExactMinovlen"]] <- usearchAttemptExactMinovlen
  #Remove trailing "/" if present
  exp_metadata[["output_path"]] <- gsub("/$", "", outputPath)
  exp_metadata[["project_name"]] <- projectName
  exp_metadata[["wildtypeSequence"]] <- wildtypeSequence
  exp_metadata[["transLibrary"]] <- transLibrary
  exp_metadata[["num_cores"]] <- numCores

  #First pipeline stage to run
  first_stage <- startStage
  last_stage <- stopStage

  ###########################
  ### MAIN
  ###########################

  ### Output file path, working and temp directories
  ###########################

  #Create working directory (if doesn't already exist)
  exp_metadata[["project_path"]] <- file.path(exp_metadata[["output_path"]], exp_metadata[["project_name"]])
  suppressWarnings(dir.create(exp_metadata[["project_path"]]))
  #Set working directory
  setwd(exp_metadata[["project_path"]])
  #Create temp directory (if doesn't already exist)
  exp_metadata[["tmp_path"]] <- file.path(exp_metadata[["project_path"]], "tmp")
  suppressWarnings(dir.create(exp_metadata[["tmp_path"]]))

  ### Get experiment design
  ###########################

  #TODO: check if all fastq files exist
  exp_metadata[["exp_design"]] <- get_experiment_design(exp_metadata)

  ### Get barcode design (if provided)
  ###########################

  if(!is.null(exp_metadata[["barcode_design_path"]])){
    exp_metadata[["barcode_design"]] <- read.table(exp_metadata[["barcode_design_path"]], header = T, stringsAsFactors = F, sep="\t")
  }

  ### Pipeline stages
  ###########################

  ### Step 0: Start pipeline tracking
  pipeline <- list()
  pipeline[['0_original']] <- exp_metadata

  ### Step 1: Run demultiplex on all fastq files
  pipeline[['1_demultiplex']] <- dimsum_stage_demultiplex(dimsum_meta = pipeline[['0_original']], demultiplex_outpath = file.path(pipeline[['0_original']][["tmp_path"]], "demultiplex"), 
    execute = (first_stage <= 1 & (last_stage == 0 | last_stage >= 1)))

  ### Step 2.1: Run FASTQC on all fastq files
  pipeline[['2_fastqc']] <- dimsum_stage_fastqc(dimsum_meta = pipeline[['1_demultiplex']], fastqc_outpath = file.path(pipeline[['1_demultiplex']][["tmp_path"]], "fastqc"), 
    execute = (first_stage <= 2 & (last_stage == 0 | last_stage >= 2)), report_outpath = file.path(pipeline[['1_demultiplex']][["project_path"]], "reports"))

  ### Step 2.2: Unzip FASTQ files if necessary
  pipeline[['2_fastq']] <- dimsum_stage_unzip(dimsum_meta = pipeline[['2_fastqc']], fastq_outpath = file.path(pipeline[['2_fastqc']][["tmp_path"]], "fastq"), 
    execute = (first_stage <= 2 & (last_stage == 0 | last_stage >= 2)))

  ### Step 2.3: Split FASTQ files
  pipeline[['2_split']] <- dimsum_stage_split(dimsum_meta = pipeline[['2_fastq']], split_outpath = file.path(pipeline[['2_fastq']][["tmp_path"]], "split"), 
    execute = (first_stage <= 2 & (last_stage == 0 | last_stage >= 2)))

  ### Step 3: Remove adapters from FASTQ files with cutadapt if necessary
  pipeline[['3_cutadapt']] <- dimsum_stage_cutadapt(dimsum_meta = pipeline[['2_split']], cutadapt_outpath = file.path(pipeline[['2_split']][["tmp_path"]], "cutadapt"), 
    execute = (first_stage <= 3 & (last_stage == 0 | last_stage >= 3)), report_outpath = file.path(pipeline[['2_split']][["project_path"]], "reports"))

  ### Step 4: Merge paired-end reads with USEARCH
  pipeline[['4_usearch']] <- dimsum_stage_usearch(dimsum_meta = pipeline[['3_cutadapt']], usearch_outpath = file.path(pipeline[['3_cutadapt']][["tmp_path"]], "usearch"), 
    execute = (first_stage <= 4 & (last_stage == 0 | last_stage >= 4)), report_outpath = file.path(pipeline[['3_cutadapt']][["project_path"]], "reports"))

  ### Step 5: Get unique aligned read counts with FASTX-Toolkit
  pipeline[['5_unique']] <- dimsum_stage_unique(dimsum_meta = pipeline[['4_usearch']], unique_outpath = file.path(pipeline[['4_usearch']][["tmp_path"]], "unique"), 
    execute = (first_stage <= 5 & (last_stage == 0 | last_stage >= 5)))

  ### Step 6: Merge variant count tables
  pipeline[['6_merge']] <- dimsum_stage_merge(dimsum_meta = pipeline[['5_unique']], merge_outpath = pipeline[['5_unique']][["project_path"]], 
    execute = (first_stage <= 6 & (last_stage == 0 | last_stage >= 6)), report_outpath = file.path(pipeline[['5_unique']][["project_path"]], "reports"))

  ### Save workspace
  ###########################

  message("\n\n\nSaving workspace image...")
  save(list = ls(environment()), file=file.path(pipeline[['6_merge']][["project_path"]], paste0(pipeline[['6_merge']][["project_name"]], '_workspace.RData')))
  message("Done")

  ### Save report html
  ###########################

  message("\n\n\nSaving summary report...")
  write(reports_summary(dimsum_meta = pipeline[['6_merge']]), file = file.path(pipeline[['6_merge']][["project_path"]], "reports_summary.html"))
  message("Done")
}


