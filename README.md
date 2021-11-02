[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/r-dimsum/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/r-dimsum/badges/version.svg)](https://anaconda.org/bioconda/r-dimsum)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/r-dimsum/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/r-dimsum)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/r-dimsum/badges/downloads.svg)](https://anaconda.org/bioconda/r-dimsum)
[![DOI](https://zenodo.org/badge/58115412.svg)](https://zenodo.org/badge/latestdoi/58115412)

<p align="left">
  <img src="./Dumpling.png" width="100">
</p>

# DiMSum

Welcome to the GitHub repository for DiMSum: An error model and pipeline for analyzing deep mutational scanning (DMS) data and diagnosing common experimental pathologies.

# Table Of Contents

1. **[Overview](#overview)**
   1. **[Pipeline Stages](docs/PIPELINE.md)**
1. **[Installation](#installation)**
1. **[Usage](#usage)**
   1. **[Command-line Arguments](docs/ARGUMENTS.md)**
   1. **[File Formats](docs/FILEFORMATS.md)**
   1. **[Demo](docs/DEMO.md)**
1. **[Bugs](#bugs)**
1. **[Citing DiMSum](#citing-dimsum)**

# Overview

The DiMSum pipeline processes raw sequencing reads (in FASTQ format) or variant counts from deep mutational scanning (DMS) experiments to calculate estimates of variant fitness (and assocated error). These estimates are suitable for use in downstream analyses of epistasis and [protein structure determination](https://github.com/lehner-lab/DMS2structure).

The DiMSum pipeline consists of five stages grouped into two modules that can be run independently:

* **_WRAP_** (DiMSum stages 1-3) processes raw FASTQ files generating a table of variant counts
* **_STEAM_** (DiMSum stages 4-5) analyses variant counts generating variant fitness and error estimates

Further details of individual DiMSum pipeline stages can be found [here](docs/PIPELINE.md).

# Installation

The easiest way to install DiMSum is by using the [bioconda package](http://bioconda.github.io/recipes/r-dimsum/README.html).
```
conda install -c bioconda r-dimsum
```

See the full [Installation Instructions](docs/INSTALLATION.md) for further details and alternative installation options.

# Usage

In the example below, DiMSum will obtain variant sequences by aligning paired-end reads in the directory "FASTQ_dir", count variant occurrences for all samples specified in the supplied [Experimental Design File](docs/FILEFORMATS.md#experimental-design-file) ("experimentDesign.txt") and calculate fitness (and error) for all variants relative to the indicated wild-type sequence.
```
DiMSum --fastqFileDir FASTQ_dir --experimentDesignPath experimentDesign.txt --wildtypeSequence AGCTAGCT
```
By default, [output files](docs/FILEFORMATS.md#output-files) are saved to the folder "DiMSum_Project" in the current working directory.

See instructions regarding [Command-line Arguments](docs/ARGUMENTS.md), [File Formats](docs/FILEFORMATS.md) and the [Demo](docs/DEMO.md) mode for full details and usage options.

# Bugs

All bug reports are highly appreciated. You may submit a bug report here on GitHub as an issue or you could send an email to ajfaure@gmail.com.

# Citing DiMSum

Please cite the following publication if you use DiMSum:

Faure, A.J., Schmiedel, J.M., Baeza-Centurion, P., Lehner B. DiMSum: an error model and pipeline for analyzing deep mutational scanning data and diagnosing common experimental pathologies. Genome Biol 21, 207 (2020). [10.1186/s13059-020-02091-3](https://doi.org/10.1186/s13059-020-02091-3)

(Vector illustration credit: <a href="https://www.vecteezy.com">Vecteezy!</a>)
