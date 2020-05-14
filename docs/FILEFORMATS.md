# Input File Formats

## Experimental design file

To run this pipeline, you will first need to describe your experimental design (e.g. in MSExcel) and save this as a tab-separated plain text file. You can download [this](./example_experimentDesign.txt) file to use as a template.

Your file must have the following columns:
* **sample_name** A sensible sample name e.g. 'input1' (alphanumeric characters only).
* **experiment_replicate** An integer identifier denoting distinct experiments (e.g. distinct plasmid library transformations) i.e. a set of input and output replicates originating from the same input biological replicate (strictly positive integer).
* **selection_id** An integer inidicating whether samples were sequenced before (0) or after (1) selection. Subsequent (serial) rounds of selection are indicated by higher numbers i.e. 2, 3, etc. (positive integer, zero included).
* **selection_replicate** An integer denoting distinct replicate selections (or biological output replicates) each derived from the same input sample (strictly positive integer). Entries should be blank (empty string) for all input samples (each input sample corresponds to a unique experiment).
* **technical_replicate** An integer denoting technical replicates (a strictly positive integer) corresponding to sample re-sequencing i.e. extracted DNA originating from the same sample split between separate sequencing lanes or files. Leave this column blank (empty string) when no technical replicates are present.
* **pair1** FASTQ file name of the first read in a given pair.
* **pair2** FASTQ file name of the second read in a given pair (omit for single-end library designs i.e. 'paired'=F).

Below is a schematic of a generic deep mutational scanning experiment indicating the corresponding entries which should be made in the experimental design file (red text). 
<p align="left">
  <img src="./DMS_experiment.png" width="600">
</p>

In addition to these mandatory columns, additional columns may be included to specify stage 2-specific options i.e. those prefixed by 'cutadapt...', which relate to constant region trimming. This allows sample-specific trimming behaviour if necessary. Options specified by columns in the experimental design file override global arguments. In the case of a growth-rate based assay, a 'generations' column can be supplied (for all output samples) in order to normalize fitness and error estimates accordingly.

## Barcode design file

If your raw FASTQ sequencing files contain multiplexed samples you will need to provide a tab-separated plain text file describing how index tags map to samples. You can download [this](./example_barcodeDesign.txt) file to use as a template.

Your file must have the following columns:
* **pair1** FASTQ file name of the first read in a given pair.
* **pair2** FASTQ file name of the second read in a given pair (omit for single-end library designs i.e. 'paired'=F).
* **barcode** Sample index tag (A/C/G/T characters only).
* **new_pair_prefix** FASTQ file prefix of demultiplexed sample reads i.e. excluding file extension (alphanumeric and underscore characters only).

When including a barcode design file, ensure that all 'new_pair_prefix' column entries correspond to 'pair1' and 'pair2' column entries in the experiment design file by appending '1.fastq' and '2.fastq' to the prefix for the first and second read respectively.

## Variant identity file

If your raw FASTQ sequencing files contain variant barcodes you will need to provide a tab-separated plain text file describing how barcodes map to variants. You can download [this](./example_variantIdentity.txt) file to use as a template. 

Your file must have the following columns:
* **barcode** DNA barcode (A/C/G/T characters only).
* **variant** Associated DNA variant (A/C/G/T characters only).
