# Overview

Welcome to the GitHub repository for DiMSum: A pipeline for pre-processing of paired-end reads from deep mutational scannning (DMS) data.

# Required Software

To run these scripts you will need the following software and associated packages:

* **[Python](https://www.python.org/downloads/) v3.6** (os, sys, argparse, biopython)
* **[R](https://www.r-project.org/) v3.3** (data.table, seqinr, optparse, parallel, reshape2, ggplot2, plyr)
* **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.3**
* **[cutadapt](https://cutadapt.readthedocs.io/en/stable/) v1.16**
* **[USEARCH 32-bit](https://drive5.com/usearch/download.html) v10.0**
* **fastx_collapser** from the [FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/download.html) v0.0.13

# Pipeline

To run this pipeline, you will first need to formalise your experiment design (e.g. in MSExcel) and save this as a tab-separated plain text file. You can download an example file to use a template [here](./example_experimentDesign.txt).

Additionally, if your raw FASTQ sequencing files contain multiplexed samples you will need to provide a tab-separated plain text file describing how barcodes map to samples. You can download an example file to use a template [here](./example_barcodeDesign.txt).

(Vector illustration credit: <a href="https://www.vecteezy.com">Vecteezy!</a>)
