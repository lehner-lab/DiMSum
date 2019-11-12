
#' dimsum__swap_reads
#'
#' Swap reads in fastq files according to adapter presence.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param cutadapt_outpath cutadapt output path (required)
#'
#' @return An updated dimsum_meta object
#' @export
dimsum__swap_reads <- function(
  dimsum_meta,
  cutadapt_outpath){
  fastq_pair_list_nunique <- dimsum_meta[['exp_design']][,c('pair1', 'pair2')]
  dimsum__swap_reads_helper <- function(
    i
    ){
    #Options for swapping read1 and read2
    temp_options_swap <- dimsum__get_cutadapt_options(dimsum_meta = dimsum_meta, exp_design_row = i, option_type = "swap")
    #Options for removing a fixed number of bases from beginning or end of either read in pair
    temp_cut_options <- dimsum__get_cutadapt_options(dimsum_meta = dimsum_meta, exp_design_row = i, option_type = "cut")
    temp_out <- system(paste0(
      "cutadapt",
      temp_options_swap,
      temp_cut_options,
      " --no-trim ",
      " --minimum-length ",
      as.character(dimsum_meta[['exp_design']][i,"cutadaptMinLength"]),
      " -e ",
      as.character(dimsum_meta[['exp_design']][i,"cutadaptErrorRate"]),
      " --untrimmed-output ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.untrimmed.fastq")),
      " --untrimmed-paired-output ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt.untrimmed.fastq")),
      " -o ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1-{name}.fastq")),
      " -p ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt1-{name}.fastq")),
      " ",
      file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[['exp_design']][i,"pair1"]),
      " ",
      file.path(dimsum_meta[["exp_design"]][i,"pair_directory"], dimsum_meta[['exp_design']][i,"pair2"]),
      " > ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1.stdout")),
      " 2> ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1.stderr"))))
    #New read1 file
    temp_out <- system(paste0(
      "cat ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt1-reverse.fastq")),
      " ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1-forward.fastq")),
      " ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt.untrimmed.fastq")),
      " > ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt2"))))
    #New read2 file
    temp_out <- system(paste0(
      "cat ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair1"], ".cutadapt1-reverse.fastq")),
      " ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt1-forward.fastq")),
      " ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt.untrimmed.fastq")),
      " > ",
      file.path(cutadapt_outpath, paste0(dimsum_meta[['exp_design']][i,"pair2"], ".cutadapt2"))))
  }
  # Setup cluster
  clust <- parallel::makeCluster(dimsum_meta[['numCores']])
  # make variables available to each core's workspace
  parallel::clusterExport(clust, list("dimsum_meta","cutadapt_outpath"), envir = environment())
  parallel::parSapply(clust,X = (1:dim(dimsum_meta[['exp_design']])[1])[!duplicated(fastq_pair_list_nunique)], dimsum__swap_reads_helper)
  parallel::stopCluster(clust)
}
