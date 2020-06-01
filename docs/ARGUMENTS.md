**[< Table Of Contents](https://github.com/lehner-lab/DiMSum#table-of-contents)**
<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# Command-line Arguments

## General

* **_--runDemo_** Run the DiMSum [Demo](DEMO.md) (default:F)
* **_--projectName_** Project name and directory where results are to be saved (default:'DiMSum_Project')
* **_--experimentDesignPath_** Path to [Experimental Design File](FILEFORMATS.md#experimental-design-file) (required, if '_--runDemo_'=F)
* **_--outputPath_** Path to directory to use for output files (default:'./' i.e. current working directory)
* **_--retainIntermediateFiles_** Should intermediate files be retained? Intermediate files can be many gigabytes, but are required to rerun DiMSum starting at intermediate pipeline stages (default:F)
* **_--startStage_** (Re-)Start DiMSum at a specific pipeline stage (default:0)
* **_--stopStage_** Stop DiMSum at a specific pipeline stage (default:5)
* **_--numCores_** Number of available CPU cores. All pipeline stages make use of parallel computing to decrease runtime if multiple cores are available (default:1)

## TRIM Arguments

* **_--cutadapt5First_** Sequence of 5' constant region to be trimmed from first (or only) read (optional)
* **_--cutadapt5Second_** Sequence of 5' constant region to be trimmed from second read in pair (optional)
* **_--cutadapt3First_** Sequence of 3' constant region to be trimmed from first (or only) read (default: reverse complement of '_--cutadapt5Second_')
* **_--cutadapt3Second_** Sequence of 3' constant region to be trimmed from second read in pair (default: reverse complement of '_--cutadapt5First_')
* **_--cutadaptMinLength_** Discard reads shorter than LENGTH after trimming (default:50)
* **_--cutadaptErrorRate_** Maximum allowed error rate for trimming constant regions (default:0.2)
* **_--cutadaptOverlap_** Minimum overlap between read and constant region for trimming (default:3)
* **_--cutadaptCut5First_** Remove fixed number of bases from start (5') of first (or only) read before constant region trimming (optional)
* **_--cutadaptCut5Second_** Remove fixed number of bases from start (5') of second read in pair before constant region trimming (optional)
* **_--cutadaptCut3First_** Remove fixed number of bases from end (3') of first (or only) read before constant region trimming (optional)
* **_--cutadaptCut3Second_** Remove fixed number of bases from end (3') of second read in pair before constant region trimming (optional)

## ALIGN Arguments

* **_--usearchMinQual_** Minimum Phred base quality score required to retain read or read pair (default:30)
* **_--usearchMaxee_** Maximum number of expected errors tolerated to retain read or read pair (default:0.5)
* **_--usearchMinlen_** Discard read (or read pair) if its length is shorter than this (default:64)
* **_--usearchMinovlen_** Discard read pair if the alignment length is shorter than this (default:16)

## PROCESS Arguments

* **_--reverseComplement_** Reverse complement variant sequences before processing? (default:F)
* **_--wildtypeSequence_** Wild-type nucleotide sequence (A/C/G/T). Lower-case bases (a/c/g/t) indicate internal constant regions to be removed (required)
* **_--permittedSequences_** Nucleotide sequence of IUPAC ambiguity codes (A/C/G/T/R/Y/S/W/K/M/B/D/H/V/N) with length matching the number of mutated positions (i.e upper-case letters) in '_--wildtypeSequence_' (default:N i.e. any substitution mutation allowed)
* **_--sequenceType_** Coding potential of sequence: either 'noncoding', 'coding' or 'auto'. If the specified wild-type nucleotide sequence ('_--wildtypeSequence_') has a valid translation without a premature STOP codon, it is assumed to be 'coding' (default:auto)
* **_--mutagenesisType_** Whether mutagenesis was performed at the nucleotide or codon/amino acid level; either 'random' or 'codon' (default:random)
* **_--maxSubstitutions_** Maximum number of nucleotide or amino acid substitutions for coding or non-coding sequences respectively (default:2)
* **_--mixedSubstitutions_** For coding sequences, are nonsynonymous variants with silent/synonymous substitutions in other codons allowed? (default:F)

## ANALYSE Arguments


