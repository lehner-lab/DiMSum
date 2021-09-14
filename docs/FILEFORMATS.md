**[< Table Of Contents](https://github.com/lehner-lab/DiMSum#table-of-contents)**
<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# File Formats

* **[Experimental Design File](#experimental-design-file)**
* **[FASTQ Files](#fastq-files)**
* **[Variant Count File](#variant-count-file)**
* **[Barcode Design File](#barcode-design-file)**
* **[Variant Identity File](#variant-identity-file)**
* **[Output Files](#output-files)**

## Experimental Design File 

**REQUIRED:** DiMSum requires a table (e.g. using Microsoft Excel) describing the experimental design that has been saved as tab-separated plain text file (see [arguments](ARGUMENTS.md#general)). You can download [this](../examples/example_experimentDesign.txt) file to use as a template.

Your file must have the following columns:
* **sample_name** A sensible sample name e.g. 'input1' (alphanumeric characters only).
* **experiment_replicate** An integer identifier denoting distinct experiments (e.g. distinct plasmid library transformations) i.e. a set of input and output replicates originating from the same input biological replicate (strictly positive integer).
* **selection_id** An integer inidicating whether samples were sequenced before (0) or after (1) selection. Subsequent (serial) rounds of selection are indicated by higher numbers i.e. 2, 3, etc. (positive integer, zero included).
* **selection_replicate** (Output samples only) An integer denoting distinct replicate selections (or biological output replicates) each derived from the same input sample (strictly positive integer). Entries should be blank (empty string) for all input samples (each input sample corresponds to a unique experiment).
* **technical_replicate** An integer denoting technical replicates (a strictly positive integer) corresponding to sample re-sequencing i.e. extracted DNA originating from the same sample split between separate sequencing lanes or files. Leave this column blank (empty string) when no technical replicates are present.
* **pair1** (_WRAP_ only) [FASTQ file](#fastq-files) name of the first read in a given pair.
* **pair2** (_WRAP_ only) [FASTQ file](#fastq-files) name of the second read in a given pair (omit for single-end library designs, see [arguments](ARGUMENTS.md#fastq-files)).

Optional columns for growth-rate based assays (download template file [here](../examples/example_experimentDesign_gr.txt)):
* **generations** (Output samples only) An estimate of the number of generations in order to normalize fitness and error estimates accordingly.
* **cell_density** An estimate of the cell density (optical density or similar) in order to estimate variant growth rates.
* **selection_time** (Output samples only) The selection time in hours in order to estimate variant growth rates.

Below is a schematic of a generic deep mutational scanning experiment indicating the corresponding entries which should be made in the experimental design file (red text). 
<p align="left">
  <img src="../DMS_experiment.png" width="600">
</p>

In addition to these mandatory columns, additional columns may be included to specify [Stage 2](PIPELINE.md#stage-2-trim-constant-regions-wrap)-specific options (see [arguments](ARGUMENTS.md#trim-arguments)), which relate to constant region trimming. This allows sample-specific trimming behaviour if necessary. Options specified by columns in the experimental design file override global arguments provided on the command-line.

## FASTQ Files

**OPTIONAL:** If processing of raw sequencing reads is required (with *WRAP*), DiMSum requires FASTQ files saved in a common directory and with a common file extension (see [arguments](ARGUMENTS.md#fastq-files)). Either FASTQ files or a [Variant Count File](#variant-count-file) can be supplied (not both).

## Variant Count File

**OPTIONAL:** If raw sequencing reads have already been processed independently of DiMSum, processing and analysis of variant counts (with DiMSum *STEAM*) requires a table (e.g. using Microsoft Excel) with variant sequences and counts for all samples (see [arguments](ARGUMENTS.md#custom-variant-count-file)). You can download [this](../examples/example_variantCounts.txt) file to use as a template. Either [FASTQ Files](#fastq-files) or a Variant Count File can be supplied (not both).

## Barcode Design File

**OPTIONAL:** If [FASTQ files](#fastq-files) contain multiplexed samples, DiMSum requires a table (e.g. using Microsoft Excel) describing how index tags map to samples that has been saved as tab-separated plain text file (see [arguments](ARGUMENTS.md#barcoded-library-design)). You can download [this](../examples/example_barcodeDesign.txt) file to use as a template.

Your file must have the following columns:
* **pair1** [FASTQ file](#fastq-files) name of the first read in a given pair.
* **pair2** [FASTQ file](#fastq-files) name of the second read in a given pair (omit for single-end library designs, see [arguments](ARGUMENTS.md#fastq-files)).
* **barcode** Sample index tag (A/C/G/T characters only).
* **new_pair_prefix** [FASTQ file](#fastq-files) prefix of demultiplexed sample reads i.e. excluding file extension (alphanumeric and underscore characters only).

When including a Barcode Design File, ensure that all 'new_pair_prefix' column entries correspond to 'pair1' and 'pair2' column entries in the [Experimental Design File](#experimental-design-file) by appending '1.fastq' and '2.fastq' to the prefix for the first and second read respectively.

## Variant Identity File

**OPTIONAL:** If the supplied sequences (supplied in the [FASTQ Files](#fastq-files) or [Variant Count File](#variant-count-file)) contain variant barcodes, DiMSum requires a table (e.g. using Microsoft Excel) describing how barcodes map to variants that has been saved as tab-separated plain text file (see [arguments](ARGUMENTS.md#barcoded-library-design)). You can download [this](../examples/example_variantIdentity.txt) file to use as a template. 

Your file must have the following columns:
* **barcode** DNA barcode (A/C/G/T characters only).
* **variant** Associated DNA variant (A/C/G/T characters only).

## Output Files

Primary output files:

* **report.html** DiMSum pipeline summary report and diagnostic plots in html format.
* **DiMSum_Project_fitness_replicates.RData** R data object with replicate (and merged) fitness scores and associated errors.
* **DiMSum_Project_variant_data_merge.RData** R data object with variant counts and statistics.

Additional output files:

* **fitness_wildtype.txt** Wild-type fitness score and associated error.
* **fitness_singles.txt** Single amino acid or nucleotide substitution variant fitness scores and associated errors.
* **fitness_doubles.txt** Double amino acid or nucleotide substitution variant fitness scores and associated errors.
* **fitness_silent.txt** Silent (synonymous) substitution variant fitness scores and associated errors (for coding sequences only).
* **fitness_singles_MaveDB.csv** [MaveDB](https://www.mavedb.org/) compatible .csv file with single amino acid or nucleotide substitution variant fitness scores and associated errors.
* **DiMSum_Project_variant_data_merge.tsv** Tab-separated plain text file with variant counts and statistics.
* **DiMSum_Project_nobarcode_variant_data_merge.tsv** Tab-separated plain text file with sequenced barcodes that were not found in the variant identity file.
* **DiMSum_Project_indel_variant_data_merge.tsv** Tab-separated plain text file with rejected indel variants.
* **DiMSum_Project_rejected_variant_data_merge.tsv** Tab-separated plain text file with remaining rejected variants (internal constant region mutants, mutations inconsistent with the library design or variants with too many substitutions).
