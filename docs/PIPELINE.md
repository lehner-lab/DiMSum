**[< Table Of Contents](https://github.com/lehner-lab/DiMSum#table-of-contents)**
<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# Pipeline Stages

* [Stage 0: **DEMULTIPLEX** raw reads](#stage-0-demultiplex-raw-reads-wrap)
* [Stage 1: **QC** raw reads](#stage-1-qc-raw-reads-wrap)
* [Stage 2: **TRIM** constant regions](#stage-2-trim-constant-regions-wrap)
* [Stage 3: **ALIGN** paired-end reads](#stage-3-align-paired-end-reads-wrap)
* [Stage 4: **PROCESS** variants](#stage-4-process-variants-steam)
* [Stage 5: **ANALYSE** counts](#stage-5-analyse-counts-steam)

## Stage 0: **DEMULTIPLEX** raw reads (_WRAP_)

Demultiplex samples and trim read barcodes using *[Cutadapt](INSTALLATION.md)* (optional). This stage is run if a [Barcode Design File](FILEFORMATS.md#barcode-design-file) is supplied (see [arguments](ARGUMENTS.md#multiplexed-fastq-files)).

## Stage 1: **QC** raw reads (_WRAP_)

Produce raw read quality reports using *[FastQC](INSTALLATION.md)* (and unzip and split FASTQ files if necessary).

## Stage 2: **TRIM** constant regions (_WRAP_)

Remove constant region sequences from read 5’ and 3’ ends using *[Cutadapt](INSTALLATION.md)*. By default the sequences of 3' constant regions are assumed to be the reverse complement of 5' constant region sequences (see [stage-specific arguments](ARGUMENTS.md#trim-arguments)).

## Stage 3: **ALIGN** paired-end reads (_WRAP_)

Align overlapping read pairs using *[VSEARCH](INSTALLATION.md)* and filter resulting variants according to base quality, expected number of errors and constituent read length (see [stage-specific arguments](ARGUMENTS.md#align-arguments)). Unique variant sequences are then tallied using *[Starcode](INSTALLATION.md)*. For [Trans library designs](ARGUMENTS.md#trans-library-design), read pairs are simply concatenated. For single-end libraries, reads are only filtered.

## Stage 4: **PROCESS** variants (_STEAM_)

Combine sample-wise variant counts and statistics to produce a unified results data.table. After aggregating counts across technical replicates, variants are processed and filtered according to user specifications (see [stage-specific arguments](ARGUMENTS.md#process-arguments)):
* **4.1** For [Barcoded library designs](ARGUMENTS.md#barcoded-library-design), read counts are aggregated at the variant level for barcode/variant mappings specified in the [Variant Identity File](FILEFORMATS.md#variant-identity-file). Undefined/misread barcodes are ignored.
* **4.2** Indel variants (defined as those not matching the wild-type nucleotide sequence length) are removed if necessary (see ['_--indels_' argument](ARGUMENTS.md#process-arguments)).
* **4.3** If internal constant region(s) are specified, these are excised from all substitution variants if a perfect match is found (see ['_--wildtypeSequence_' argument](ARGUMENTS.md#process-arguments)).
* **4.4** Substitution variants with mutations inconsistent with the library design are removed (see ['_--permittedSequences_' argument](ARGUMENTS.md#process-arguments)).
* **4.5** Substitution variants with more substitutions than desired are also removed (see ['_--maxSubstitutions_' argument](ARGUMENTS.md#process-arguments)).
* **4.6** Finally, nonsynonymous substitution variants with synonymous substitutions in other codons are removed if necessary (see ['_--mixedSubstitutions_' argument](ARGUMENTS.md#process-arguments)).

## Stage 5: **ANALYSE** counts (_STEAM_)

Calculate fitness and error estimates for a user-specified subset of variants (see [stage-specific arguments](ARGUMENTS.md#analyse-arguments)):
* **5.1** Optionally remove low count variants according to user-specified soft/hard thresholds to minimise the impact of "fictional" variants from sequencing errors.
* **5.2** Calculate replicate normalisation parameters (scale and shift) to minimise inter-replicate fitness differences.
* **5.3** Fit the error model to a high confidence subset of variants to determine additive and multiplicative error terms.
* **5.4** Aggregate variant fitness and error at the amino acid level if the target molecule is a coding sequence.
* **5.5** Optionally normalise fitness and error estimates by the number of generations in the case of a growth-rate based assay (see [Experiment Design File](FILEFORMATS.md#experimental-design-file)).
* **5.6** Merge fitness scores between replicates in a weighted manner that takes into account their respective errors.
