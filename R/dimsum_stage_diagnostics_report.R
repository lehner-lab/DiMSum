
#' dimsum_stage_diagnostics_report
#'
#' Generate diagnostic plots for all samples.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param variant_data path to variant data R Data file (required)
#' @param report_outpath diagnostics report output path (required)
#' @param input_samples_pattern input samples string pattern (default: "^input")
#' @param output_samples_pattern diagnostics report output path (default: "^output")
#'
#' @return Nothing
#' @export
dimsum_stage_diagnostics_report <- function(
  dimsum_meta,
  variant_data,
  report_outpath,
  input_samples_pattern = "^input",
  output_samples_pattern = "^output"
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))
  #load variant data from RData file
  load(variant_data)
  #Exclude indel variants
  vdm_noindel <- variant_data_merge[Nins_nt==0 & Ndel_nt==0,]
  #Sample names
  input_samples <- colnames(vdm_noindel)[grep(input_samples_pattern, colnames(vdm_noindel))]
  output_samples <- colnames(vdm_noindel)[grep(output_samples_pattern, colnames(vdm_noindel))]
  #Plot 1: bimodality of sequencing counts in input (per # nucleotide mutations)
  #Histogram of input counts split by number of nucleotide mutations (restricting to 1-4)
  if(length(input_samples)!=0){
    dimsum__sample_count_distributions_by_ntmut(input_dt = vdm_noindel[data.table::between(Nmut_nt,1,4),.SD,.SDcols=c(input_samples, "Nmut_nt", "Nmut_aa")],
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_count_hist_input.png"),
      title = paste0("Substitution variant count distributions (input samples)"), sequence_type = dimsum_meta[["sequenceType"]])
    dimsum__sample_count_distributions_by_ntmut(input_dt = vdm_noindel[data.table::between(Nmut_nt,1,4),.SD,.SDcols=c(input_samples, "Nmut_nt", "Nmut_aa")],
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_count_hist_input.pdf"),
      title = paste0("Substitution variant count distributions (input samples)"), sequence_type = dimsum_meta[["sequenceType"]])
  }
  #Histogram of output counts split by number of nucleotide mutations (restricting to 1-4)
  if(length(output_samples)!=0){
    dimsum__sample_count_distributions_by_ntmut(input_dt = vdm_noindel[data.table::between(Nmut_nt,1,4),.SD,.SDcols=c(output_samples, "Nmut_nt", "Nmut_aa")],
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_count_hist_output.png"),
      title = paste0("Substitution variant count distributions (output samples)"), sequence_type = dimsum_meta[["sequenceType"]])
  }
  #Plot 2: 'flaps' in replicate versus replicate plots (indicating some variants are 'highly' present in one replicate but just background noise in another)
  #Plot pairwise input sample count correlations (restricting to 1-2 nucleotide mutations)
  if(length(input_samples)!=0){
    temp_dt <- vdm_noindel[,.SD,.SDcols = c(input_samples, "Nmut_nt")]
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nmut_nt,1,2),.SD,.SDcols = input_samples]+1), 
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_scatterplotmatrix_input.png"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      title = "Substitution variant inter-sample count correlations (input samples)")
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nmut_nt,1,2),.SD,.SDcols = input_samples]+1), 
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_scatterplotmatrix_input.pdf"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      title = "Substitution variant inter-sample count correlations (input samples)")
  }
  #Plot pairwise output sample count correlations
  if(length(output_samples)!=0){
    temp_dt <- vdm_noindel[,.SD,.SDcols = c(output_samples, "Nmut_nt")]
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nmut_nt,1,2),.SD,.SDcols = output_samples]+1), 
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_scatterplotmatrix_output.png"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      title = "Substitution variant inter-sample count correlations (output samples)")
  }
  #Plot pairwise output sample count correlations
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- vdm_noindel[,.SD,.SDcols = c(input_samples, output_samples, "Nmut_nt")]
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nmut_nt,1,2),.SD,.SDcols = c(input_samples, output_samples)]+1), 
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_scatterplotmatrix_all.png"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      title = "Substitution variant inter-sample count correlations (all samples)",
      cut = as.factor(temp_dt[data.table::between(Nmut_nt,1,2),Nmut_nt]),
      size = 0.1)
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nmut_nt,1,2),.SD,.SDcols = c(input_samples, output_samples)]+1), 
      output_file = file.path(report_outpath, "dimsum_stage_diagnostics_report_scatterplotmatrix_all.pdf"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      title = "Substitution variant inter-sample count correlations (all samples)",
      cut = as.factor(temp_dt[data.table::between(Nmut_nt,1,2),Nmut_nt]),
      size = 0.1)
  }
}

