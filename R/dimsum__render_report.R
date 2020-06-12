
#' dimsum__render_report
#'
#' Generate a summary report in html format.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param initialise whether to initialise plots with touch (default:F)
#' @param finalise whether to finalse the report with end time (default:F)
#'
#' @return Nothing
#' @export
dimsum__render_report <- function(
  dimsum_meta,
  initialise=FALSE,
  finalise=FALSE
  ){

  #Initialise plots
  if(initialise){
    all_files <- list(
      "dimsum__fastqc_report_pair1_fastqc.png",
      "dimsum__fastqc_report_pair2_fastqc.png",
      "dimsum__cutadapt_report_pair1.png",
      "dimsum__cutadapt_report_pair2.png",
      "dimsum__vsearch_report_paircounts.png",
      "dimsum__vsearch_report_mergedlength.png",
      "dimsum__merge_report_nucmutationpercentages.png",
      "dimsum__merge_report_nucmutationcounts.png",
      "dimsum__merge_report_aamutationpercentages.png",
      "dimsum__merge_report_aamutationcounts.png",
      "dimsum__diagnostics_report_count_hist_input_nt.png",
      "dimsum__diagnostics_report_count_hist_input_aa.png",
      "dimsum__diagnostics_report_scatterplotmatrix_all.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.png",
      "dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.png",
      "dimsum_stage_fitness_report_1_errormodel_repspec.png",
      "dimsum_stage_fitness_report_1_errormodel_leaveoneout_qqplot.png")
    lapply(all_files, function(x){write(NULL, file = file.path(dimsum_meta[["project_path"]], "reports", x))})
  }

  #Set end time
  if(finalise){
    #Set run end time
    dimsum_meta[["end_time"]] <- Sys.time()
  }

  #Document settings
  doc_settings <- list()
  doc_settings[["dimsum_version"]] <- as.character(packageVersion("DiMSum"))
  doc_settings[["project_name"]] <- dimsum_meta[["projectName"]]
  doc_settings[["start_time"]] <- dimsum_meta[["start_time"]]
  doc_settings[["end_time"]] <- dimsum_meta[["end_time"]]
  doc_settings[["arg_list"]] <- dimsum_meta[["arg_list"]]
  doc_settings[["show_qualitycontrol1"]] <- TRUE
  doc_settings[["show_qualitycontrol2"]] <- TRUE
  doc_settings[["show_trim1"]] <- TRUE
  doc_settings[["show_trim2"]] <- TRUE
  doc_settings[["show_align1"]] <- TRUE
  doc_settings[["show_process1"]] <- TRUE
  doc_settings[["show_processA"]] <- TRUE
  doc_settings[["show_processC"]] <- TRUE
  doc_settings[["show_analyse1"]] <- TRUE
  doc_settings[["show_analyseN"]] <- TRUE

  #Remove unnecessary items
  #Noncoding sequence
  if(dimsum_meta[["sequenceType"]]!="coding"){
    doc_settings[["show_processA"]] <- FALSE
    doc_settings[["show_processC"]] <- FALSE
  }
  #Random mutagenesis
  if(dimsum_meta[["mutagenesisType"]]=="random"){
    doc_settings[["show_processC"]] <- FALSE
  }
  #Single-end sequencing
  if(!dimsum_meta[["paired"]]){
    doc_settings[["show_qualitycontrol2"]] <- FALSE
    doc_settings[["show_trim2"]] <- FALSE
  }
  #No fitness normalisation
  if(!dimsum_meta[["fitnessNormalise"]]){
    doc_settings[["show_analyseN"]] <- FALSE
  }
  #No error model or No biological replicates (requirement for normalisation and error model fit)
  if(length(unique(dimsum_meta[["exp_design"]][,"experiment_replicate"]))<2 | !dimsum_meta[["fitnessErrorModel"]]){
    doc_settings[["show_analyse1"]] <- FALSE
    doc_settings[["show_analyseN"]] <- FALSE
  }
  #STEAM reports only
  if(!is.null(dimsum_meta[["countPath"]])){
    doc_settings[["show_qualitycontrol1"]] <- FALSE
    doc_settings[["show_qualitycontrol2"]] <- FALSE
    doc_settings[["show_trim1"]] <- FALSE
    doc_settings[["show_trim2"]] <- FALSE
    doc_settings[["show_align1"]] <- FALSE
  }

  #Save document settings
  save(doc_settings, file = file.path(dimsum_meta[["project_path"]], "reports", "report_settings.RData"))

  #Copy report R markdown
  file.copy(
    from = system.file("rmd", "report.Rmd", package = "DiMSum"),
    to = file.path(dimsum_meta[["project_path"]], "reports", "report.Rmd"))

  #Copy DiMSum image
  file.copy(
    from = system.file("rmd", "Dumpling.png", package = "DiMSum"),
    to = file.path(dimsum_meta[["project_path"]], "reports", "Dumpling.png"))

  #Render R markdown
  markdown_file <- file.path(dimsum_meta[["project_path"]], "reports", "report.Rmd")
  { sink("/dev/null")
    suppressMessages(rmarkdown::render(
      input = markdown_file,
      output_file = "report.html",
      output_dir = dimsum_meta[["project_path"]]))
    sink()
  }

}
