#!/usr/bin/env Rscript

#Software version
version_info <- "DiMSum_v0.1"
#Display software version
message(version_info)

###########################
### CHECK DEPENDENCIES
###########################

#Binaries
required_binaries <- c(
  "cutadapt", 
  "cat", 
  "cp", 
  "fastqc", 
  "head", 
  "fastx_collapser", 
  "gunzip", 
  "usearch")
which_binaries <- Sys.which(required_binaries)
missing_binaries <- names(which_binaries)[which_binaries==""]
if(length(missing_binaries)!=0){
  stop(paste0("Required executables not installed. Please install the following software: ", paste(missing_binaries, sep = ", ")), call. = FALSE)
}

#R packages
required_packages <- c(
  "data.table", 
  "seqinr", 
  "ShortRead", 
  "parallel", 
  "reshape2", 
  "ggplot2", 
  "plyr", 
  "GGally",
  "hexbin",
  "optparse")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)!=0){
  stop(paste0("Required R packages not installed. Please install the following packages: ", paste(missing_packages, sep = ", ")), call. = FALSE)
}

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
  make_option(opt_str=c("--cutadaptCut5First"), type="integer", help = "cutadapt: remove bases from start of first read (before constant region trimming)"),
  make_option(opt_str=c("--cutadaptCut5Second"), type="integer", help = "cutadapt: remove bases from start of second read (before constant region trimming)"),
  make_option(opt_str=c("--cutadaptCut3First"), type="integer", help = "cutadapt: remove bases from end of first read (before constant region trimming)"),
  make_option(opt_str=c("--cutadaptCut3Second"), type="integer", help = "cutadapt: remove bases from end of second read (before constant region trimming)"),
  make_option(opt_str=c("--cutadapt5First"), help = "cutadapt: sequence of 5' constant region (of the first read)"),
  make_option(opt_str=c("--cutadapt5Second"), help = "cutadapt: sequence of 5' constant region (of the second read)"),
  make_option(opt_str=c("--cutadapt3First"), help = "cutadapt: sequence of 3' constant region (of the first read; default: reverse complement of cutadapt5Second)"),
  make_option(opt_str=c("--cutadapt3Second"), help = "cutadapt: sequence of 3' constant region (of the second read; default: reverse complement of cutadapt5First)"),
  make_option(opt_str=c("--cutadaptMinLength", "-n"), type="integer", default=50, help = "cutadapt: Discard reads shorter than LENGTH (default:50)"),
  make_option(opt_str=c("--cutadaptErrorRate", "-a"), type="double", default=0.2, help = "cutadapt: Maximum allowed error rate (default:0.2)"),
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
### LOAD PACKAGES
###########################

suppressWarnings(suppressMessages(require(data.table)))
suppressWarnings(suppressMessages(library(seqinr)))
suppressWarnings(suppressMessages(library(ShortRead)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(GGally)))
suppressWarnings(suppressMessages(library(hexbin)))

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
### RUN
###########################

dimsum(
  fastqFileDir=arg_list[["fastqFileDir"]],
  fastqFileExtension=arg_list[["fastqFileExtension"]],
  gzipped=arg_list[["gzipped"]],
  stranded=arg_list[["stranded"]],
  barcodeDesignPath=arg_list[["barcodeDesignPath"]],
  barcodeErrorRate=arg_list[["barcodeErrorRate"]],
  experimentDesignPath=arg_list[["experimentDesignPath"]],
  cutadaptCut5First=arg_list[["cutadaptCut5First"]],
  cutadaptCut5Second=arg_list[["cutadaptCut5Second"]],
  cutadaptCut3First=arg_list[["cutadaptCut3First"]],
  cutadaptCut3Second=arg_list[["cutadaptCut3Second"]],
  cutadapt5First=arg_list[["cutadapt5First"]],
  cutadapt5Second=arg_list[["cutadapt5Second"]],
  cutadapt3First=arg_list[["cutadapt3First"]],
  cutadapt3Second=arg_list[["cutadapt3Second"]],
  cutadaptMinLength=arg_list[["cutadaptMinLength"]],
  cutadaptErrorRate=arg_list[["cutadaptErrorRate"]],
  usearchMinQual=arg_list[["usearchMinQual"]],
  usearchMaxee=arg_list[["usearchMaxee"]],
  usearchMinlen=arg_list[["usearchMinlen"]],
  usearchMinovlen=arg_list[["usearchMinovlen"]],
  usearchAttemptExactMinovlen=arg_list[["usearchAttemptExactMinovlen"]],
  outputPath=arg_list[["outputPath"]],
  projectName=arg_list[["projectName"]],
  wildtypeSequence=arg_list[["wildtypeSequence"]],
  transLibrary=arg_list[["transLibrary"]],
  startStage=arg_list[["startStage"]],
  stopStage=arg_list[["stopStage"]],
  numCores=arg_list[["numCores"]])
