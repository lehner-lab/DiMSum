**[< Table Of Contents](https://github.com/lehner-lab/DiMSum#table-of-contents)**
<p align="left">
  <img src="../Dumpling.png" width="100">
</p>

# Command-line Arguments

## General

* **_--runDemo_** Run the DiMSum [Demo](DEMO.md) (default:F)
* **_--projectName_** Project name and directory where results are to be saved (default:'DiMSum_Project')
* **_--experimentDesignPath_** Path to [Experimental Design File](FILEFORMATS.md#experimental-design-file) (required, if '_--runDemo_'=F)
* **_--wildtypeSequence_** Wild-type nucleotide sequence (A/C/G/T). Lower-case letters (a/c/g/t) indicate internal constant regions to be removed during _STEAM_ (see [PROCESS Arguments](#process-arguments)).
* **_--outputPath_** Path to directory to use for output files (default:'./' i.e. current working directory)
* **_--retainIntermediateFiles_** Should intermediate files be retained? Intermediate files can be many gigabytes, but are required to rerun DiMSum starting at intermediate pipeline stages (default:F)
* **_--startStage_** (Re-)Start DiMSum at a specific pipeline stage (default:0)
* **_--stopStage_** Stop DiMSum at a specific pipeline stage (default:5)
* **_--numCores_** Number of available CPU cores. All pipeline stages make use of parallel computing to decrease runtime if multiple cores are available (default:1)

## TRIM Arguments



## ALIGN Arguments


## PROCESS Arguments


## ANALYSE Arguments


