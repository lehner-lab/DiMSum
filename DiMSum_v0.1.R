#!/usr/bin/env Rscript

###########################
### CHECK DEPENDENCIES
###########################

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
  optparse::make_option(opt_str=c("--paired"), type="logical", default=T, help = "Is the library design paired-end?"),
  optparse::make_option(opt_str=c("--barcodeDesignPath", "-b"), help = "Path to barcode design file (tab-separated plain text file with barcode design)"),
  optparse::make_option(opt_str=c("--barcodeErrorRate"), type="double", default=0.25, help = "Maximum allowed error rate for the barcode"),
  optparse::make_option(opt_str=c("--experimentDesignPath", "-e"), help = "Path to experimental design file (tab-separated plain text file with replicate structure)"),
  optparse::make_option(opt_str=c("--experimentDesignPairDuplicates"), type="logical", default=F, help = "Are duplicate FASTQ files permitted in experimental design file?"),
  optparse::make_option(opt_str=c("--cutadaptCut5First"), type="integer", help = "cutadapt: remove bases from start of first read (before constant region trimming)"),
  optparse::make_option(opt_str=c("--cutadaptCut5Second"), type="integer", help = "cutadapt: remove bases from start of second read (before constant region trimming)"),
  optparse::make_option(opt_str=c("--cutadaptCut3First"), type="integer", help = "cutadapt: remove bases from end of first read (before constant region trimming)"),
  optparse::make_option(opt_str=c("--cutadaptCut3Second"), type="integer", help = "cutadapt: remove bases from end of second read (before constant region trimming)"),
  optparse::make_option(opt_str=c("--cutadapt5First"), help = "cutadapt: sequence of 5' constant region (of the first read)"),
  optparse::make_option(opt_str=c("--cutadapt5Second"), help = "cutadapt: sequence of 5' constant region (of the second read)"),
  optparse::make_option(opt_str=c("--cutadapt3First"), help = "cutadapt: sequence of 3' constant region (of the first read; default: reverse complement of cutadapt5Second)"),
  optparse::make_option(opt_str=c("--cutadapt3Second"), help = "cutadapt: sequence of 3' constant region (of the second read; default: reverse complement of cutadapt5First)"),
  optparse::make_option(opt_str=c("--cutadaptMinLength", "-n"), type="integer", default=50, help = "cutadapt: Discard reads shorter than LENGTH"),
  optparse::make_option(opt_str=c("--cutadaptErrorRate", "-a"), type="double", default=0.2, help = "cutadapt: Maximum allowed error rate"),
  optparse::make_option(opt_str=c("--usearchMinQual", "-q"), type="integer", help = "USEARCH: minimum observed base quality to retain read pair"),
  optparse::make_option(opt_str=c("--usearchMaxee", "-m"), type="double", help = "USEARCH: maximum number of expected errors to retain read pair"),
  optparse::make_option(opt_str=c("--usearchMinlen"), type="integer", default=64, help = "USEARCH: Discard pair if either read is shorter than this"),
  optparse::make_option(opt_str=c("--usearchMinovlen"), type="integer", default=16, help = "USEARCH: discard pair if alignment is shorter than given value"),
  optparse::make_option(opt_str=c("--outputPath", "-o"), help = "Path to directory to use for output files"),
  optparse::make_option(opt_str=c("--projectName", "-p"), help = "Project name"),
  optparse::make_option(opt_str=c("--wildtypeSequence", "-w"), help = "Wild-type nucleotide sequence (A/C/G/T). Lower-case letters (a/c/g/t) indicate internal constant regions to be removed before fitness calculations"),
  optparse::make_option(opt_str=c("--permittedSequences"), help = "A sequence of nucleotide codes (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N) with length matching the number of mutated positions i.e upper-case letters in wild-type nucleotide sequence (default:any base at mutated positions)"),
  optparse::make_option(opt_str=c("--sequenceType", "-u"), default="auto", help = "Coding potential of sequence; either noncoding/coding/auto (default:auto)"),
  optparse::make_option(opt_str=c("--transLibrary", "-r"), type="logical", default=F, help = "Trans library design i.e. read pairs correspond to distinct peptides (no overlap)"),
  optparse::make_option(opt_str=c("--bayesianDoubleFitness", "-y"), type="logical", default=F, help = "Fitness: improve double mutant fitness estimates using Bayesian framework (DISABLED: still in development)"),
  optparse::make_option(opt_str=c("--bayesianDoubleFitnessLamD", "-d"), type="double", default=0.025, help = "Fitness: Poisson distribution for score likelihood (default:0.025)"),
  optparse::make_option(opt_str=c("--fitnessMinInputCountAll"), type="integer", default=0, help = "Fitness: minimum input read count (in all replicate) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessMinInputCountAny"), type="integer", default=0, help = "Fitness: minimum input read count (in any replicate) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessHighConfidenceCount"), type="integer", default=10, help = "Fitness: minimum mean input read count for high confidence variants (default:10)"),
  optparse::make_option(opt_str=c("--fitnessNormalise"), type="logical", default=T, help = "Normalise fitness values to minimise inter-replicate differences (default:T)"),
  optparse::make_option(opt_str=c("--fitnessDoubleHighConfidenceCount"), type="integer", default=50, help = "Fitness: minimum input replicate read count for doubles used to derive prior for Bayesian doubles correction (default:50)"),
  optparse::make_option(opt_str=c("--maxSubstitutions"), type="integer", default=2, help = "Maximum number of nucleotide or amino acid substitutions for coding or non-coding sequences respectively (default:2)"),
  optparse::make_option(opt_str=c("--retainIntermediateFiles", "-f"), type="logical", default=F, help = "Should intermediate files be retained? (default:F)"),
  optparse::make_option(opt_str=c("--splitChunkSize", "-k"), type="double", default=3758096384, help = "FASTQ file split chunk size in bytes (default:3758096384)"),
  optparse::make_option(opt_str=c("--retainedReplicates"), default='all', help = "Comma-separated list of Input replicates (or experiment ids) to retain (default:'all')"),
  optparse::make_option(opt_str=c("--startStage", "-s"), type="integer", default=1, help = "Start at a specified pipeline stage"),
  optparse::make_option(opt_str=c("--stopStage", "-t"), type="integer", default=0, help = "Stop at a specified pipeline stage (default: 0 i.e. no stop condition)"),
  optparse::make_option(opt_str=c("--numCores", "-c"), type="integer", default=1, help = "Number of available CPU cores")  
)

arg_list <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
#Display arguments
message(paste("\n\n\n*******", "DiMSum wrapper command-line arguments", "*******\n\n\n"))
message(paste(formatDL(unlist(arg_list)), collapse = "\n"))

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
  paired=arg_list[["paired"]],
  barcodeDesignPath=arg_list[["barcodeDesignPath"]],
  barcodeErrorRate=arg_list[["barcodeErrorRate"]],
  experimentDesignPath=arg_list[["experimentDesignPath"]],
  experimentDesignPairDuplicates=arg_list[["experimentDesignPairDuplicates"]],
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
  outputPath=arg_list[["outputPath"]],
  projectName=arg_list[["projectName"]],
  wildtypeSequence=arg_list[["wildtypeSequence"]],
  permittedSequences=arg_list[["permittedSequences"]],
  sequenceType=arg_list[["sequenceType"]],
  transLibrary=arg_list[["transLibrary"]],
  bayesianDoubleFitness=arg_list[["bayesianDoubleFitness"]],
  bayesianDoubleFitnessLamD=arg_list[["bayesianDoubleFitnessLamD"]],
  fitnessMinInputCountAll=arg_list[["fitnessMinInputCountAll"]],
  fitnessMinInputCountAny=arg_list[["fitnessMinInputCountAny"]],
  fitnessHighConfidenceCount=arg_list[["fitnessHighConfidenceCount"]],
  fitnessDoubleHighConfidenceCount=arg_list[["fitnessDoubleHighConfidenceCount"]],
  maxSubstitutions=arg_list[["maxSubstitutions"]],
  retainIntermediateFiles=arg_list[["retainIntermediateFiles"]],
  splitChunkSize=arg_list[["splitChunkSize"]],
  retainedReplicates=arg_list[["retainedReplicates"]],
  startStage=arg_list[["startStage"]],
  stopStage=arg_list[["stopStage"]],
  numCores=arg_list[["numCores"]])
