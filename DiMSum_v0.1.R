#!/usr/bin/env Rscript

#Software version
version_info <- "DiMSum_v0.1"
#Display software version
message(version_info)

###########################
### COMMAND-LINE OPTIONS
###########################

suppressWarnings(suppressMessages(library(optparse)))

option_list <- list(
  make_option(opt_str=c("--fastqFileDir", "-i"), help = "Path to directory with input FASTQ files"),
  make_option(opt_str=c("--fastqFileExtension", "-l"), default='.fastq', help = "FASTQ file extension"),
  make_option(opt_str=c("--gzipped", "-g"), type="logical", default=T, help = "Are FASTQ files are gzipped?"),
  make_option(opt_str=c("--stranded"), type="logical", default=T, help = "Is the library design stranded?"),
  make_option(opt_str=c("--barcodeDesignPath", "-b"), help = "Path to barcode design file (tab-separated plain text file with barcode design)"),
  make_option(opt_str=c("--barcodeErrorRate"), type="double", default=0.25, help = "Maximum allowed error rate for the barcode (default:0.25)"),
  make_option(opt_str=c("--experimentDesignPath", "-e"), help = "Path to experimental design file (tab-separated plain text file with replicate structure)"),
  make_option(opt_str=c("--cutadaptCut5First"), type="integer", help = "cutadapt: remove bases from start of first read (before adapter trimming)"),
  make_option(opt_str=c("--cutadaptCut5Second"), type="integer", help = "cutadapt: remove bases from start of second read (before adapter trimming)"),
  make_option(opt_str=c("--cutadaptCut3First"), type="integer", help = "cutadapt: remove bases from end of first read (before adapter trimming)"),
  make_option(opt_str=c("--cutadaptCut3Second"), type="integer", help = "cutadapt: remove bases from end of second read (before adapter trimming)"),
  make_option(opt_str=c("--cutadapt5First"), help = "cutadapt: sequence of an adapter ligated to the 5' end (of the first read)"),
  make_option(opt_str=c("--cutadapt5Second"), help = "cutadapt: sequence of an adapter ligated to the 5' end (of the second read)"),
  make_option(opt_str=c("--cutadapt3First"), help = "cutadapt: sequence of an adapter ligated to the 3' end (of the first read)"),
  make_option(opt_str=c("--cutadapt3Second"), help = "cutadapt: sequence of an adapter ligated to the 3' end (of the second read)"),
  make_option(opt_str=c("--cutadaptMinLength", "-n"), type="integer", default=50, help = "cutadapt: Discard reads shorter than LENGTH (default:50)"),
  make_option(opt_str=c("--cutadaptErrorRate", "-a"), type="double", default=0.2, help = "cutadapt: Maximum allowed error rate (default:0.2)"),
  make_option(opt_str=c("--cutadaptDiscardUntrimmed"), type="logical", default=F, help = "cutadapt: Discard untrimmed read pairs (default:F)"),
  make_option(opt_str=c("--usearchMinQual", "-q"), type="integer", help = "USEARCH: minimum observed base quality to retain read pair"),
  make_option(opt_str=c("--usearchMaxee", "-m"), type="double", help = "USEARCH: maximum number of expected errors to retain read pair"),
  make_option(opt_str=c("--usearchMinlen"), type="integer", default=64, help = "USEARCH: Discard pair if either read is shorter than this (default:64)"),
  make_option(opt_str=c("--usearchMinovlen"), type="integer", default=16, help = "USEARCH: discard pair if alignment is shorter than given value (default:16)"),
  make_option(opt_str=c("--usearchAttemptExactMinovlen"), type="logical", default=F, help = "USEARCH: Attempt exact alignment of --usearchMinovlen (default:F)"),
  make_option(opt_str=c("--outputPath", "-o"), help = "Path to directory to use for output files"),
  make_option(opt_str=c("--projectName", "-p"), help = "Project name"),
  make_option(opt_str=c("--wildtypeSequence", "-w"), help = "Wild-type nucleotide sequence"),
  make_option(opt_str=c("--transLibrary", "-r"), type="logical", default=F, help = "Trans library design i.e. read pairs correspond to distinct peptides (no overlap)"),
  make_option(opt_str=c("--startStage", "-s"), type="integer", default=1, help = "Start at a specified pipeline stage"),
  make_option(opt_str=c("--stopStage", "-t"), type="integer", default=0, help = "Stop at specified pipeline stage (default:0, no stop condition)"),
  make_option(opt_str=c("--numCores", "-c"), type="integer", default=1, help = "Number of available CPU cores")  
)

arg_list <- parse_args(OptionParser(option_list=option_list))
#Display arguments
message(paste("\n\n\n*******", "DiMSum command-line arguments", "*******\n\n\n"))
print(arg_list)

###########################
### PACKAGES
###########################

suppressWarnings(suppressMessages(require(data.table)))
suppressWarnings(suppressMessages(library(seqinr)))
suppressWarnings(suppressMessages(library(ShortRead)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(GGally)))

###########################
### SCRIPTS
###########################

#Get DiMSum base directory
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
dimsum_base_dir <- normalizePath(dirname(script.name))

#Source scripts
filelist = list.files(file.path(dimsum_base_dir, 'scripts/'))
x <- sapply(paste0(file.path(dimsum_base_dir, 'scripts/'),filelist),source,.GlobalEnv)

###########################
### GLOBALS
###########################

#Metadata object
exp_metadata <- list()

### Save metadata
#Remove trailing "/" if present
exp_metadata[["fastq_path_original"]] <- gsub("/$", "", arg_list$fastqFileDir)
exp_metadata[["fastq_file_extension"]] <- arg_list$fastqFileExtension
exp_metadata[["gzipped"]] <- arg_list$gzipped
exp_metadata[["stranded"]] <- arg_list$stranded
exp_metadata[["barcode_design_path"]] <- arg_list$barcodeDesignPath
exp_metadata[["barcodeErrorRate"]] <- arg_list$barcodeErrorRate
exp_metadata[["experiment_design_path"]] <- arg_list$experimentDesignPath
exp_metadata[["cutadaptCut5First"]] <- arg_list$cutadaptCut5First
exp_metadata[["cutadaptCut5Second"]] <- arg_list$cutadaptCut5Second
exp_metadata[["cutadaptCut3First"]] <- arg_list$cutadaptCut3First
exp_metadata[["cutadaptCut3Second"]] <- arg_list$cutadaptCut3Second
exp_metadata[["cutadapt5First"]] <- arg_list$cutadapt5First
exp_metadata[["cutadapt5Second"]] <- arg_list$cutadapt5Second
exp_metadata[["cutadapt3First"]] <- arg_list$cutadapt3First
exp_metadata[["cutadapt3Second"]] <- arg_list$cutadapt3Second
exp_metadata[["cutadaptMinLength"]] <- arg_list$cutadaptMinLength
exp_metadata[["cutadaptErrorRate"]] <- arg_list$cutadaptErrorRate
exp_metadata[["cutadaptDiscardUntrimmed"]] <- arg_list$cutadaptDiscardUntrimmed
exp_metadata[["usearchMinQual"]] <- arg_list$usearchMinQual
exp_metadata[["usearchMaxee"]] <- arg_list$usearchMaxee
exp_metadata[["usearchMinlen"]] <- arg_list$usearchMinlen
exp_metadata[["usearchMinovlen"]] <- arg_list$usearchMinovlen
exp_metadata[["usearchAttemptExactMinovlen"]] <- arg_list$usearchAttemptExactMinovlen
#Remove trailing "/" if present
exp_metadata[["output_path"]] <- gsub("/$", "", arg_list$outputPath)
exp_metadata[["project_name"]] <- arg_list$projectName
exp_metadata[["wildtypeSequence"]] <- arg_list$wildtypeSequence
exp_metadata[["transLibrary"]] <- arg_list$transLibrary
exp_metadata[["num_cores"]] <- arg_list$numCores

#First pipeline stage to run
first_stage <- arg_list$startStage
last_stage <- arg_list$stopStage

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

#TODO: check if file exists
#TODO: check if all fastq files exist
exp_metadata[["exp_design"]] <- read.table(exp_metadata[["experiment_design_path"]], header = T, stringsAsFactors = F, sep="\t")
#Add original FASTQ directory
exp_metadata[["exp_design"]]$pair_directory <- exp_metadata[["fastq_path_original"]]

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

### Step 2: Run FASTQC on all fastq files
pipeline[['2_fastqc']] <- dimsum_stage_fastqc(dimsum_meta = pipeline[['1_demultiplex']], fastqc_outpath = file.path(pipeline[['1_demultiplex']][["tmp_path"]], "fastqc"), 
  execute = (first_stage <= 2 & (last_stage == 0 | last_stage >= 2)), report_outpath = file.path(pipeline[['1_demultiplex']][["project_path"]], "reports"))

### Step 3: Unzip FASTQ files if necessary
pipeline[['3_fastq']] <- dimsum_stage_unzip(dimsum_meta = pipeline[['2_fastqc']], fastq_outpath = file.path(pipeline[['2_fastqc']][["tmp_path"]], "fastq"), 
  execute = (first_stage <= 3 & (last_stage == 0 | last_stage >= 3)))

### Step 4: Split FASTQ files
pipeline[['4_split']] <- dimsum_stage_split(dimsum_meta = pipeline[['3_fastq']], split_outpath = file.path(pipeline[['3_fastq']][["tmp_path"]], "split"), 
  execute = (first_stage <= 4 & (last_stage == 0 | last_stage >= 4)))

### Step 5: Remove adapters from FASTQ files with cutadapt if necessary
pipeline[['5_cutadapt']] <- dimsum_stage_cutadapt(dimsum_meta = pipeline[['4_split']], cutadapt_outpath = file.path(pipeline[['4_split']][["tmp_path"]], "cutadapt"), 
  execute = (first_stage <= 5 & (last_stage == 0 | last_stage >= 5)), report_outpath = file.path(pipeline[['4_split']][["project_path"]], "reports"))

### Step 6: Merge paired-end reads with USEARCH
pipeline[['6_usearch']] <- dimsum_stage_usearch(dimsum_meta = pipeline[['5_cutadapt']], usearch_outpath = file.path(pipeline[['5_cutadapt']][["tmp_path"]], "usearch"), 
  execute = (first_stage <= 6 & (last_stage == 0 | last_stage >= 6)), report_outpath = file.path(pipeline[['5_cutadapt']][["project_path"]], "reports"))

### Step 7: Get unique aligned read counts with FASTX-Toolkit
pipeline[['7_unique']] <- dimsum_stage_unique(dimsum_meta = pipeline[['6_usearch']], unique_outpath = file.path(pipeline[['6_usearch']][["tmp_path"]], "unique"), 
  execute = (first_stage <= 7 & (last_stage == 0 | last_stage >= 7)))

### Step 8: Merge variant count tables
pipeline[['8_merge']] <- dimsum_stage_merge(dimsum_meta = pipeline[['7_unique']], merge_outpath = pipeline[['7_unique']][["project_path"]], 
  execute = (first_stage <= 8 & (last_stage == 0 | last_stage >= 8)), report_outpath = file.path(pipeline[['7_unique']][["project_path"]], "reports"))

### Save workspace
###########################

message("\n\n\nSaving workspace image...")
save.image(file=file.path(pipeline[['8_merge']][["project_path"]], paste0(pipeline[['8_merge']][["project_name"]], '_workspace.RData')))
message("Done")

### Save report html
###########################

message("\n\n\nSaving summary report...")
write(gsub("PROJECT_NAME", pipeline[['8_merge']][["project_name"]], reports_summary_template), file = file.path(pipeline[['8_merge']][["project_path"]], "reports_summary.html"))
message("Done")


