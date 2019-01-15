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
  "cat", 
  "cp", 
  "cutadapt", 
  "fastqc", 
  "fastx_collapser", 
  "gunzip", 
  "head", 
  "usearch")
which_binaries <- Sys.which(required_binaries)
missing_binaries <- names(which_binaries)[which_binaries==""]
if(length(missing_binaries)!=0){
  stop(paste0("Required executables not installed. Please install the following software: ", paste(missing_binaries, sep = ", ")), call. = FALSE)
}

#R packages
required_packages <- c(
  "DiMSum")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)!=0){
  stop(paste0("Required R packages not installed. Please install the following packages: ", paste(missing_packages, sep = ", ")), call. = FALSE)
}

###########################
### COMMAND-LINE OPTIONS
###########################

option_list <- list(
  optparse::make_option(opt_str=c("--fastqFileDir", "-i"), help = "Path to directory with input FASTQ files"),
  optparse::make_option(opt_str=c("--fastqFileExtension", "-l"), default='.fastq', help = "FASTQ file extension"),
  optparse::make_option(opt_str=c("--gzipped", "-g"), type="logical", default=T, help = "Are FASTQ files are gzipped?"),
  optparse::make_option(opt_str=c("--stranded"), type="logical", default=T, help = "Is the library design stranded?"),
  optparse::make_option(opt_str=c("--barcodeDesignPath", "-b"), help = "Path to barcode design file (tab-separated plain text file with barcode design)"),
  optparse::make_option(opt_str=c("--barcodeErrorRate"), type="double", default=0.25, help = "Maximum allowed error rate for the barcode (default:0.25)"),
  optparse::make_option(opt_str=c("--experimentDesignPath", "-e"), help = "Path to experimental design file (tab-separated plain text file with replicate structure)"),
  optparse::make_option(opt_str=c("--cutadaptCut5First"), type="integer", help = "cutadapt: remove bases from start of first read (before constant region trimming)"),
  optparse::make_option(opt_str=c("--cutadaptCut5Second"), type="integer", help = "cutadapt: remove bases from start of second read (before constant region trimming)"),
  optparse::make_option(opt_str=c("--cutadaptCut3First"), type="integer", help = "cutadapt: remove bases from end of first read (before constant region trimming)"),
  optparse::make_option(opt_str=c("--cutadaptCut3Second"), type="integer", help = "cutadapt: remove bases from end of second read (before constant region trimming)"),
  optparse::make_option(opt_str=c("--cutadapt5First"), help = "cutadapt: sequence of 5' constant region (of the first read)"),
  optparse::make_option(opt_str=c("--cutadapt5Second"), help = "cutadapt: sequence of 5' constant region (of the second read)"),
  optparse::make_option(opt_str=c("--cutadapt3First"), help = "cutadapt: sequence of 3' constant region (of the first read; default: reverse complement of cutadapt5Second)"),
  optparse::make_option(opt_str=c("--cutadapt3Second"), help = "cutadapt: sequence of 3' constant region (of the second read; default: reverse complement of cutadapt5First)"),
  optparse::make_option(opt_str=c("--cutadaptMinLength", "-n"), type="integer", default=50, help = "cutadapt: Discard reads shorter than LENGTH (default:50)"),
  optparse::make_option(opt_str=c("--cutadaptErrorRate", "-a"), type="double", default=0.2, help = "cutadapt: Maximum allowed error rate (default:0.2)"),
  optparse::make_option(opt_str=c("--usearchMinQual", "-q"), type="integer", help = "USEARCH: minimum observed base quality to retain read pair"),
  optparse::make_option(opt_str=c("--usearchMaxee", "-m"), type="double", help = "USEARCH: maximum number of expected errors to retain read pair"),
  optparse::make_option(opt_str=c("--usearchMinlen"), type="integer", default=64, help = "USEARCH: Discard pair if either read is shorter than this (default:64)"),
  optparse::make_option(opt_str=c("--usearchMinovlen"), type="integer", default=16, help = "USEARCH: discard pair if alignment is shorter than given value (default:16)"),
  optparse::make_option(opt_str=c("--usearchAttemptExactMinovlen"), type="logical", default=F, help = "USEARCH: Attempt exact alignment of --usearchMinovlen (default:F)"),
  optparse::make_option(opt_str=c("--outputPath", "-o"), help = "Path to directory to use for output files"),
  optparse::make_option(opt_str=c("--projectName", "-p"), help = "Project name"),
  optparse::make_option(opt_str=c("--wildtypeSequence", "-w"), help = "Wild-type nucleotide sequence"),
  optparse::make_option(opt_str=c("--transLibrary", "-r"), type="logical", default=F, help = "Trans library design i.e. read pairs correspond to distinct peptides (no overlap)"),
  optparse::make_option(opt_str=c("--startStage", "-s"), type="integer", default=1, help = "Start at a specified pipeline stage (default:1)"),
  optparse::make_option(opt_str=c("--stopStage", "-t"), type="integer", default=0, help = "Stop at a specified pipeline stage (default:0 i.e. no stop condition)"),
  optparse::make_option(opt_str=c("--numCores", "-c"), type="integer", default=1, help = "Number of available CPU cores")  
)

arg_list <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
#Display arguments
message(paste("\n\n\n*******", "DiMSum command-line arguments", "*******\n\n\n"))
print(arg_list)

###########################
### LOAD DIMSUM PACKAGE
###########################

library(DiMSum)

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
