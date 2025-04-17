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
### FUNCTIONS
###########################

status_message <- function(...){cat(sprintf(...), sep='', file=stdout())}

###########################
### COMMAND-LINE OPTIONS
###########################

option_list <- list(
  optparse::make_option(opt_str=c("--runDemo"), type="logical", default=F, help = "Run the DiMSum demo (default:F)"),
  optparse::make_option(opt_str=c("--fastqFileDir", "-i"), default='./', help = "Path to directory containing input FASTQ files (default:'./' i.e. current working directory)"),
  optparse::make_option(opt_str=c("--fastqFileExtension", "-l"), default='.fastq', help = "FASTQ file extension (default:'.fastq')"),
  optparse::make_option(opt_str=c("--gzipped", "-g"), type="logical", default=T, help = "Are FASTQ files are gzipped? (default:T)"),
  optparse::make_option(opt_str=c("--stranded"), type="logical", default=T, help = "Is the library design stranded? (default:T)"),
  optparse::make_option(opt_str=c("--paired"), type="logical", default=T, help = "Is the library design paired-end? (default:T)"),
  optparse::make_option(opt_str=c("--barcodeDesignPath", "-b"), help = "Path to barcode design file (tab-separated plain text file with barcode design)"),
  optparse::make_option(opt_str=c("--barcodeErrorRate"), type="double", default=0.25, help = "Maximum allowed error rate for barcode to be matched (default:0.25)"),
  optparse::make_option(opt_str=c("--experimentDesignPath", "-e"), help = "Path to Experimental Design File (required if '--runDemo'=F)"),
  optparse::make_option(opt_str=c("--experimentDesignPairDuplicates"), type="logical", default=F, help = "Are multiple instances of FASTQ files in the Experimental Design File permitted? (default:F)"),
  optparse::make_option(opt_str=c("--barcodeIdentityPath"), help = "Path to Variant Identity File (tab-separated plain text file mapping barcodes to variants)"),
  optparse::make_option(opt_str=c("--countPath"), help = "Path to Variant Count File for analysis with STEAM only (tab-separated plain text file with sample counts for all variants)"),
  optparse::make_option(opt_str=c("--synonymSequencePath"), help = "Path to Synonym Sequences File with coding sequences for which synonymous variant fitness should be quantified (default: plain text file with one coding nucleotide sequence per line)"),
  optparse::make_option(opt_str=c("--cutadaptCut5First"), type="integer", help = "Remove fixed number of bases from start (5') of first (or only) read before constant region trimming (optional)"),
  optparse::make_option(opt_str=c("--cutadaptCut5Second"), type="integer", help = "Remove fixed number of bases from start (5') of second read in pair before constant region trimming (optional)"),
  optparse::make_option(opt_str=c("--cutadaptCut3First"), type="integer", help = "Remove fixed number of bases from end (3') of first (or only) read before constant region trimming (optional)"),
  optparse::make_option(opt_str=c("--cutadaptCut3Second"), type="integer", help = "Remove fixed number of bases from end (3') of second read in pair before constant region trimming (optional)"),
  optparse::make_option(opt_str=c("--cutadapt5First"), help = "Sequence of 5' constant region to be trimmed from first (or only) read (optional)"),
  optparse::make_option(opt_str=c("--cutadapt5Second"), help = "Sequence of 5' constant region to be trimmed from second read in pair (optional)"),
  optparse::make_option(opt_str=c("--cutadapt3First"), help = "Sequence of 3' constant region to be trimmed from first (or only) read (default: reverse complement of '--cutadapt5Second')"),
  optparse::make_option(opt_str=c("--cutadapt3Second"), help = "Sequence of 3' constant region to be trimmed from second read in pair (default: reverse complement of '--cutadapt5First')"),
  optparse::make_option(opt_str=c("--cutadaptMinLength", "-n"), type="integer", default=50, help = "Discard reads shorter than LENGTH after trimming (default:50)"),
  optparse::make_option(opt_str=c("--cutadaptErrorRate", "-a"), type="double", default=0.2, help = "Maximum allowed error rate for trimming constant regions (default:0.2)"),
  optparse::make_option(opt_str=c("--cutadaptOverlap"), type="integer", default=3, help = "Minimum overlap between read and constant region for trimming (default:3)"),
  optparse::make_option(opt_str=c("--vsearchMinQual", "-q"), type="integer", default=30, help = "Minimum Phred base quality score required to retain read or read pair (default:30)"),
  optparse::make_option(opt_str=c("--vsearchMaxQual"), type="integer", default=41, help = "Maximum Phred base quality score accepted when reading (and used when writing) FASTQ files; cannot be greater than 93 (default:41)"),
  optparse::make_option(opt_str=c("--vsearchMaxee", "-m"), type="double", default=0.5, help = "Maximum number of expected errors tolerated to retain read or read pair (default:0.5)"),
  optparse::make_option(opt_str=c("--vsearchMinovlen"), type="integer", default=10, help = "Discard read pair if the alignment length is shorter than this (default:10)"),
  optparse::make_option(opt_str=c("--outputPath", "-o"), default='./', help = "Path to directory to use for output files (default:'./' i.e. current working directory)"),
  optparse::make_option(opt_str=c("--projectName", "-p"), default='DiMSum_Project', help = "Project name and directory where results are to be saved (default:'DiMSum_Project')"),
  optparse::make_option(opt_str=c("--wildtypeSequence", "-w"), help = "Wild-type nucleotide sequence (A/C/G/T). Lower-case bases (a/c/g/t) indicate internal constant regions to be removed (required if '--runDemo'=F)"),
  optparse::make_option(opt_str=c("--permittedSequences"), help = "Nucleotide sequence of IUPAC ambiguity codes (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N) with length matching the number of mutated positions (i.e upper-case letters) in '--wildtypeSequence' (default:N i.e. any substitution mutation allowed)"),
  optparse::make_option(opt_str=c("--reverseComplement"), type="logical", default=F, help = "Reverse complement sequence (default:F)"),
  optparse::make_option(opt_str=c("--sequenceType", "-u"), default="auto", help = "Coding potential of sequence: either 'noncoding', 'coding' or 'auto'. If the specified wild-type nucleotide sequence ('--wildtypeSequence') has a valid translation without a premature STOP codon, it is assumed to be 'coding' (default:'auto')"),
  optparse::make_option(opt_str=c("--mutagenesisType"), default="random", help = "Whether mutagenesis was performed at the nucleotide or codon/amino acid level; either 'random' or 'codon' (default:'random')"),
  optparse::make_option(opt_str=c("--transLibrary", "-r"), type="logical", default=F, help = "Paired-end reads correspond to distinct molecules? (default:F)"),
  optparse::make_option(opt_str=c("--transLibraryReverseComplement"), type="logical", default=F, help = "Reverse complement second read in pair (default:F)"),
  optparse::make_option(opt_str=c("--bayesianDoubleFitness", "-y"), type="logical", default=F, help = "In development: improve double mutant fitness estimates using Bayesian framework (DISABLED: still in development)"),
  optparse::make_option(opt_str=c("--bayesianDoubleFitnessLamD", "-d"), type="double", default=0.025, help = "In development: Poisson distribution for score likelihood (default:0.025)"),
  optparse::make_option(opt_str=c("--fitnessMinInputCountAll"), default="0", help = "Minimum input read count (in all replicates) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessMinInputCountAny"), default="0", help = "Minimum input read count (in any replicate) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessMinOutputCountAll"), default="0", help = "Minimum output read count (in all replicates) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessMinOutputCountAny"), default="0", help = "Minimum output read count (in any replicates) to be retained during fitness calculations (default:0)"),
  optparse::make_option(opt_str=c("--fitnessHighConfidenceCount"), type="integer", default=10, help = "In development: minimum mean input read count for high confidence variants (default:10)"),
  optparse::make_option(opt_str=c("--fitnessDoubleHighConfidenceCount"), type="integer", default=50, help = "In development: minimum input replicate read count for doubles used to derive prior for Bayesian doubles correction (default:50)"),
  optparse::make_option(opt_str=c("--fitnessNormalise"), type="logical", default=T, help = "Normalise fitness values to minimise inter-replicate differences (default:T)"),
  optparse::make_option(opt_str=c("--fitnessErrorModel"), type="logical", default=T, help = "Fit fitness error model (default:T)"),
  optparse::make_option(opt_str=c("--fitnessDropoutPseudocount"), type="integer", default=0, help = "Pseudocount added to output replicates with dropout i.e. variants present in input but absent from output (default:0)"),
  optparse::make_option(opt_str=c("--indels"), default="none", help = "Indel variants to be retained: either 'all', 'none' or a comma-separated list of sequence lengths (default:'none')"),
  optparse::make_option(opt_str=c("--maxSubstitutions"), type="integer", default=2, help = "Maximum number of nucleotide or amino acid substitutions for coding or non-coding sequences respectively (default:2)"),
  optparse::make_option(opt_str=c("--mixedSubstitutions"), type="logical", default=F, help = "For coding sequences, are nonsynonymous variants with silent/synonymous substitutions in other codons allowed? (default:F)"),
  optparse::make_option(opt_str=c("--retainIntermediateFiles", "-f"), type="logical", default=F, help = "Should intermediate files be retained? Intermediate files can be many gigabytes, but are required to rerun DiMSum starting at intermediate pipeline stages (default:F)"),
  optparse::make_option(opt_str=c("--splitChunkSize", "-k"), type="double", default=3758096384, help = "Internal: FASTQ file split chunk size in bytes (default:3758096384)"),
  optparse::make_option(opt_str=c("--retainedReplicates"), default='all', help = "Comma-separated list of (integer) experiment replicates to retain or 'all' (default:'all')"),
  optparse::make_option(opt_str=c("--startStage", "-s"), type="integer", default=0, help = "(Re-)Start DiMSum at a specific pipeline stage (default:0)"),
  optparse::make_option(opt_str=c("--stopStage", "-t"), type="integer", default=5, help = "Stop DiMSum at a specific pipeline stage (default:5)"),
  optparse::make_option(opt_str=c("--numCores", "-c"), type="integer", default=1, help = "Number of available CPU cores. All pipeline stages make use of parallel computing to decrease runtime if multiple cores are available (default:1)")  
)

arg_list <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
#Display arguments
status_message(paste("\n\n\n*******", "DiMSum wrapper command-line arguments", "*******\n\n\n"))
status_message(paste(formatDL(unlist(arg_list)), collapse = "\n"))

###########################
### LOAD DIMSUM PACKAGE
###########################

library(DiMSum)

###########################
### RUN
###########################

dimsum(
  runDemo=arg_list[["runDemo"]],
  fastqFileDir=arg_list[["fastqFileDir"]],
  fastqFileExtension=arg_list[["fastqFileExtension"]],
  gzipped=arg_list[["gzipped"]],
  stranded=arg_list[["stranded"]],
  paired=arg_list[["paired"]],
  barcodeDesignPath=arg_list[["barcodeDesignPath"]],
  barcodeErrorRate=arg_list[["barcodeErrorRate"]],
  experimentDesignPath=arg_list[["experimentDesignPath"]],
  experimentDesignPairDuplicates=arg_list[["experimentDesignPairDuplicates"]],
  barcodeIdentityPath=arg_list[["barcodeIdentityPath"]],
  countPath=arg_list[["countPath"]],
  synonymSequencePath=arg_list[["synonymSequencePath"]],
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
  cutadaptOverlap=arg_list[["cutadaptOverlap"]],
  vsearchMinQual=arg_list[["vsearchMinQual"]],
  vsearchMaxQual=arg_list[["vsearchMaxQual"]],
  vsearchMaxee=arg_list[["vsearchMaxee"]],
  vsearchMinovlen=arg_list[["vsearchMinovlen"]],
  outputPath=arg_list[["outputPath"]],
  projectName=arg_list[["projectName"]],
  wildtypeSequence=arg_list[["wildtypeSequence"]],
  permittedSequences=arg_list[["permittedSequences"]],
  reverseComplement=arg_list[["reverseComplement"]],
  sequenceType=arg_list[["sequenceType"]],
  mutagenesisType=arg_list[["mutagenesisType"]],
  transLibrary=arg_list[["transLibrary"]],
  transLibraryReverseComplement=arg_list[["transLibraryReverseComplement"]],
  bayesianDoubleFitness=arg_list[["bayesianDoubleFitness"]],
  bayesianDoubleFitnessLamD=arg_list[["bayesianDoubleFitnessLamD"]],
  fitnessMinInputCountAll=arg_list[["fitnessMinInputCountAll"]],
  fitnessMinInputCountAny=arg_list[["fitnessMinInputCountAny"]],
  fitnessMinOutputCountAll=arg_list[["fitnessMinOutputCountAll"]],
  fitnessMinOutputCountAny=arg_list[["fitnessMinOutputCountAny"]],
  fitnessHighConfidenceCount=arg_list[["fitnessHighConfidenceCount"]],
  fitnessDoubleHighConfidenceCount=arg_list[["fitnessDoubleHighConfidenceCount"]],
  fitnessNormalise=arg_list[["fitnessNormalise"]],
  fitnessErrorModel=arg_list[["fitnessErrorModel"]],
  fitnessDropoutPseudocount=arg_list[["fitnessDropoutPseudocount"]],
  indels=arg_list[["indels"]],
  maxSubstitutions=arg_list[["maxSubstitutions"]],
  mixedSubstitutions=arg_list[["mixedSubstitutions"]],
  retainIntermediateFiles=arg_list[["retainIntermediateFiles"]],
  splitChunkSize=arg_list[["splitChunkSize"]],
  retainedReplicates=arg_list[["retainedReplicates"]],
  startStage=arg_list[["startStage"]],
  stopStage=arg_list[["stopStage"]],
  numCores=arg_list[["numCores"]])
