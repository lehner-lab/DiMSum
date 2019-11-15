
#' dimsum_stage_diagnostics_report
#'
#' Generate diagnostic plots for all samples.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param report_outpath diagnostics report output path (required)
#' @param input_samples_pattern input samples string pattern (default: "_s0_")
#' @param output_samples_pattern diagnostics report output path (default: "_s1_")
#'
#' @return Nothing
#' @export
dimsum_stage_diagnostics_report <- function(
  dimsum_meta,
  report_outpath,
  input_samples_pattern = "_s0_",
  output_samples_pattern = "_s1_"
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))
  #load variant data from RData file
  load(dimsum_meta[["variant_data_merge_path"]])

  #Sample names
  input_samples <- colnames(variant_data_merge)[grep(input_samples_pattern, colnames(variant_data_merge))]
  output_samples <- colnames(variant_data_merge)[grep(output_samples_pattern, colnames(variant_data_merge))]

  #Plot 1: bimodality of sequencing counts in input (per # nucleotide mutations)
  #Histogram of input counts split by number of nucleotide mutations (restricting to 1-4)
  if(length(input_samples)!=0){
    dimsum__sample_count_distributions_by_ntham(input_dt = variant_data_merge[data.table::between(Nham_nt,1,4),.SD,.SDcols=c(input_samples, "Nham_nt", "Nham_aa")],
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_count_hist_input.png"),
      title = paste0("Substitution variant count distributions (input samples)"), sequence_type = dimsum_meta[["sequenceType"]])
    dimsum__sample_count_distributions_by_ntham(input_dt = variant_data_merge[data.table::between(Nham_nt,1,4),.SD,.SDcols=c(input_samples, "Nham_nt", "Nham_aa")],
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_count_hist_input.pdf"),
      title = paste0("Substitution variant count distributions (input samples)"), sequence_type = dimsum_meta[["sequenceType"]])
  }
  #Histogram of output counts split by number of nucleotide mutations (restricting to 1-4)
  if(length(output_samples)!=0){
    dimsum__sample_count_distributions_by_ntham(input_dt = variant_data_merge[data.table::between(Nham_nt,1,4),.SD,.SDcols=c(output_samples, "Nham_nt", "Nham_aa")],
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_count_hist_output.png"),
      title = paste0("Substitution variant count distributions (output samples)"), sequence_type = dimsum_meta[["sequenceType"]])
  }
  
  #Plot 2: all-vs-all sample count correlations
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- variant_data_merge[,.SD,.SDcols = c(input_samples, output_samples, "Nham_nt")]
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nham_nt,1,2),.SD,.SDcols = c(input_samples, output_samples)]+1), 
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_scatterplotmatrix_all.png"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      title = "Single and double substitution variant inter-sample count correlations (all samples)",
      cut = as.factor(temp_dt[data.table::between(Nham_nt,1,2),Nham_nt]),
      size = 0.1,
      thresholds = list(
        "dotted" = log10(dimsum_meta[["fitnessMinInputCountAny"]] + 1), 
        "dashed" = log10(dimsum_meta[["fitnessMinInputCountAll"]] + 1)))
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nham_nt,1,2),.SD,.SDcols = c(input_samples, output_samples)]+1), 
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_scatterplotmatrix_all.pdf"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      title = "Single and double substitution variant inter-sample count correlations (all samples)",
      cut = as.factor(temp_dt[data.table::between(Nham_nt,1,2),Nham_nt]),
      size = 0.1,
      thresholds = list(
        "dotted" = log10(dimsum_meta[["fitnessMinInputCountAny"]] + 1), 
        "dashed" = log10(dimsum_meta[["fitnessMinInputCountAll"]] + 1)))
  }
}

