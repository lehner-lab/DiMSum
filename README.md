<p align="left">
  <img src="./Dumpling.png" width="100">
</p>

# Overview

Welcome to the GitHub repository for DiMSum: A pipeline for pre-processing deep mutational scanning (DMS) data.

# Required Software

To run the DiMSum pipeline you will need the following software and associated packages:

* **[R](https://www.r-project.org/) >=v3.5.2** (data.table, ggplot2, gridExtra, hexbin, optparse, parallel, plyr, reshape2, seqinr, GGally, ShortRead)
* **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.3**
* **[cutadapt](https://cutadapt.readthedocs.io/en/stable/) v2.4**
* **[USEARCH 32-bit](https://drive5.com/usearch/download.html) v10.0**
* **[starcode](https://github.com/gui11aume/starcode) v1.3**

# Installation and loading

Open R and enter:

```
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("lehner-lab/DiMSum")

# Load
library(DiMSum)

# Help
?dimsum
```

# DiMSum command-line tool

Clone the DiMSum repository and install the R package locally. The * must be replaced by what is actually downloaded and built.

```
git clone https://github.com/lehner-lab/DiMSum.git
R CMD build DiMSum
R CMD INSTALL DiMSum_*.tar.gz
```
Add the cloned DiMSum repository base directory to your path. You can do this by adding the following line at the bottom of your ~/.bashrc file:
```
export PATH=CLONED_DIMSUM_REPOSITORY:$PATH
```
Get a description of DiMSum command-line arguments with the following:
```
DiMSum -h
```

# Pipeline

The DiMSum pipeline processes raw sequencing reads (in FASTQ format) from deep mutational scanning (DMS) experiments to produce variant counts, fitness and error estimates. These estimates are suitable for use in downstream analyses of epistasis and [protein structure determination](https://github.com/lehner-lab/DMS2structure).

## Stage 1: DEMULTIPLEX READS

Demultiplex samples and trim read barcodes using cutadapt (optional). This stage is run if a barcode design file is supplied (see 'barcodeDesignPath' argument). Stage-specific arguments: 'barcodeDesignPath' and 'barcodeErrorRate'.

## Stage 2: ASSESS READ QUALITY

Produce raw read quality reports using FastQC.

## Stage 3: TRIM CONSTANT REGIONS

Remove constant region sequences from read 5’ and 3’ ends using cutadapt. 5' constant region sequences are required (3' sequences are optional). By default the sequences of 3' constant regions are assumed to be the reverse complement of 5' constant region sequences. Stage-specific arguments: 'cutadaptCut5First', 'cutadaptCut5Second', 'cutadaptCut3First', 'cutadaptCut3Second', 'cutadapt5First', 'cutadapt5Second', 'cutadapt3First', 'cutadapt3Second', 'cutadaptMinLength', 'cutadaptErrorRate'.

## Stage 4: ALIGN PAIRED-END READS

Align overlapping read pairs using USEARCH (paired-end cis libraries only i.e. 'paired'=T, 'transLibrary'=F) or alternatively concatenate read pairs (paired-end trans libraries only i.e. 'transLibrary'=T), and filter resulting variants according to base quality, expected number of errors and constituent read length (including those from single-end libraries i.e. 'paired'=F). Stage-specific arguments: 'usearchMinQual', 'usearchMaxee', 'usearchMinlen', 'usearchMinovlen'.

## Stage 5: COUNT UNIQUE VARIANTS

Tally counts of unique variant sequences using [starcode](https://github.com/gui11aume/starcode).

## Stage 6: MERGE SAMPLE STATISTICS

Combine sample-wise variant counts and statistics to produce a unified results data.table. Variant counts are aggregated across technical replicates.

## Stage 7: CALCULATE FITNESS

Calculate fitness and error estimates for a user-specified subset of substitution variants:
* **7.1** Indel variants not matching the WT sequence length are discarded.
* **7.2** If internal constant region(s) are specified, these are  removed from all variants if a perfect match is found.
* **7.3** An error model is fit to all variants to determine the contributions of count-based (Poisson), replicate and over-sequencing error terms.
* **7.4** Nucleotide variants with low input read counts likely representing sequencing errors ('fitnessMinInputCountAll', 'fitnessMinInputCountAny') are filtered out. Additionally, variants exceeding the maximum number of nucleotide or amino acid substitutions for coding or non-coding sequences respectively ('fitnessMaxSubstitutions') are filtered out. For amino acid substitutions, nonsynonymous variants with synonymous variants in other codons are discarded.
* **7.5** Variants are aggregated at the amino acid level if the target molecule is a protein ('sequenceType'=coding). 
* **7.6** Fitness and estimates of the associated error are then calculated with respect to the corresponding wild-type sequence score using the model derived in **7.3** above.

An optional step (*coming soon: still in development*) will improve double mutant fitness estimates for low frequency variants using a Bayesian approach that incorporates priors based on observed single mutant counts ('bayesianDoubleFitness', 'bayesianDoubleFitnessLamD', 'fitnessHighConfidenceCount', 'fitnessDoubleHighConfidenceCount'). In the case of a growth-rate based assay, a 'generations' column can be supplied in the experimental design file in order to normalize fitness and error estimates accordingly (see below). Finally, fitness scores are merged between replicates in a weighted manner taking into account their respective errors.

# Experimental design file

To run this pipeline, you will first need to describe your experimental design (e.g. in MSExcel) and save this as a tab-separated plain text file. You can download [this](./example_experimentDesign.txt) file to use as a template.

Your file must have the following columns:
* **sample_name** A sensible sample name e.g. 'input1' (alphanumeric characters only).
* **transformation_replicate** An integer identifier denoting distinct variant library transformations i.e. a set of input and output replicates originating from the same plasmid library transformation or input biological replicate (strictly positive integer).
* **selection_id** An integer inidicating whether samples were sequenced before (0) or after (1) selection. Subsequent (serial) rounds of selection are indicated by higher numbers i.e. 2, 3, etc. (positive integer, zero included).
* **selection_replicate** An integer denoting distinct replicate selections (or biological output replicates) each derived from the same input sample (strictly positive integer). Entries should be blank (empty string) for all input samples (each input sample corresponds to a unique transformation).
* **technical_replicate** An integer denoting technical replicates (a strictly positive integer) corresponding to sample re-sequencing i.e. extracted DNA originating from the same sample split between separate sequencing lanes or files. Leave this column blank (empty string) when no technical replicates are present.
* **pair1** FASTQ file name of the first read in a given pair.
* **pair2** FASTQ file name of the second read in a given pair (omit for single-end library designs i.e. 'paired'=F).

Below is a schematic of a generic deep mutational scanning experiment indicating the corresponding entries which should be made in the experimental design file (red text). 
<p align="left">
  <img src="./DMS_experiment.png" width="600">
</p>

In addition to these mandatory columns, additional columns may be included to specify stage 3-specific options i.e. those prefixed by 'cutadapt...', which relate to constant region trimming. This allows sample-specific trimming behaviour if necessary. Options specified by columns in the experimental design file override global arguments. In the case of a growth-rate based assay, a 'generations' column can be supplied (for all output samples) in order to normalize fitness and error estimates accordingly.

# Barcode design file

If your raw FASTQ sequencing files contain multiplexed samples you will need to provide a tab-separated plain text file describing how barcodes map to samples. You can download [this](./example_barcodeDesign.txt) file to use as a template.

Your file must have the following columns:
* **pair1** FASTQ file name of the first read in a given pair.
* **pair2** FASTQ file name of the second read in a given pair (omit for single-end library designs i.e. 'paired'=F).
* **barcode** Sample DNA barcode (A/C/G/T characters only).
* **new_pair_prefix** FASTQ file prefix of demultiplexed sample reads i.e. excluding file extension (alphanumeric and underscore characters only).

When including a barcode design file, ensure that all 'new_pair_prefix' column entries correspond to 'pair1' and 'pair2' column entries in the experiment design file by appending '1.fastq' and '2.fastq' to the prefix for the first and second read respectively.

# Output

* **PROJECT_DIR/PROJECT_NAME_variant_data_merge.RData** R data object with variant counts and statistics ('variant_data_merge' data.table).
* **PROJECT_DIR/reports_summary.html** DiMSum pipeline summary report and diagnostic plots in html format.




(Vector illustration credit: <a href="https://www.vecteezy.com">Vecteezy!</a>)
