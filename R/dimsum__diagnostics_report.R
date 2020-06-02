
#' dimsum__diagnostics_report
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
#' @import data.table
dimsum__diagnostics_report <- function(
  dimsum_meta,
  report_outpath,
  input_samples_pattern = "_s0_",
  output_samples_pattern = "_s1_"
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))

  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = dimsum_meta[["variant_data_merge_path"]],
    execute = TRUE,
    exit = FALSE)

  #load variant data from RData file
  load(dimsum_meta[["variant_data_merge_path"]])

  #Sample names
  input_samples <- colnames(variant_data_merge)[grep(input_samples_pattern, colnames(variant_data_merge))]
  output_samples <- colnames(variant_data_merge)[grep(output_samples_pattern, colnames(variant_data_merge))]

  #Max number of nucleotide substitutions to plot
  max_nsubs <- dimsum_meta[["maxSubstitutions"]]
  if(dimsum_meta[["sequenceType"]]=="coding"){
    max_nsubs <- max_nsubs*3
  }
  #Set maximum limit of 12 substitutions to display
  max_nsubs <- min(c(max_nsubs, 12))

  #Error rate
  error_rate <- 10^min(-dimsum_meta[["usearchMinQual"]]/10, log10(dimsum_meta[["usearchMaxee"]] / nchar(dimsum_meta[["wildtypeSequence"]])))

  #Plot 1: Histogram of input counts split by number of nucleotide mutations
  if(length(input_samples)!=0){
    dimsum__sample_count_distributions(
      input_dt = variant_data_merge[data.table::between(Nham_nt,0,max_nsubs),.SD,.SDcols=c(input_samples, "Nham_nt", "WT")],
      output_file_prefix = file.path(report_outpath, "dimsum__diagnostics_report_count_hist_input_nt"),
      #title = paste0("Nucleotide substitution variant count distributions (input samples)"),
      error_rate = error_rate,
      seq_length = nchar(dimsum_meta[["wildtypeSequence"]]))
  }

  #Plot 2: Histogram of input counts split by number of nucleotide mutations and number of amino acid mutations
  max_nsubs <- min(dimsum_meta[["maxSubstitutions"]], 2)
  if(dimsum_meta[["sequenceType"]]=="coding" & max_nsubs<=2){
    max_nsubs <- max_nsubs*3
    if(length(input_samples)!=0){
      dimsum__sample_count_distributions(
        input_dt = variant_data_merge[data.table::between(Nham_nt,0,max_nsubs),.SD,.SDcols=c(input_samples, "Nham_nt", "Nham_aa", "WT")],
        output_file_prefix = file.path(report_outpath, "dimsum__diagnostics_report_count_hist_input_aa"),
        #title = paste0("Nucleotide substitution variant count distributions (input samples)"),
        error_rate = error_rate,
        seq_length = nchar(dimsum_meta[["wildtypeSequence"]]), height = 12)
    }
  }

  #Plot 3a: all-vs-all sample count correlations - all variants
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- variant_data_merge[,.SD,.SDcols = c(input_samples, output_samples, "Nham_nt")]
    names(temp_dt)[grep("_count", names(temp_dt))] <- dimsum__plot_samplename(names(temp_dt)[grep("_count", names(temp_dt))])
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[,.SD,.SDcols = dimsum__plot_samplename(c(input_samples, output_samples))]+1), 
      output_file_prefix = file.path(report_outpath, "dimsum__diagnostics_report_scatterplotmatrix_all"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      #title = "Substitution variant inter-sample count correlations (all samples)",
      size = 0.1,
      thresholds = list(
        "dotted" = log10(dimsum_meta[["fitnessMinInputCountAny"]] + 1), 
        "dashed" = log10(dimsum_meta[["fitnessMinInputCountAll"]] + 1)))
  }

  #Plot 3b: all-vs-all sample count correlations - only singles and doubles
  if(length(input_samples)!=0 | length(output_samples)!=0){
    temp_dt <- variant_data_merge[,.SD,.SDcols = c(input_samples, output_samples, "Nham_nt")]
    names(temp_dt)[grep("_count", names(temp_dt))] <- dimsum__plot_samplename(names(temp_dt)[grep("_count", names(temp_dt))])
    dimsum__ggpairs_binhex(
      input_dt = log10(temp_dt[data.table::between(Nham_nt,1,2),.SD,.SDcols = dimsum__plot_samplename(c(input_samples, output_samples))]+1), 
      output_file_prefix = file.path(report_outpath, "dimsum__diagnostics_report_scatterplotmatrix_singles_doubles"),
      xlab = "log10(variant count + 1)",
      ylab = "log10(variant count + 1)",
      #title = "Single and double substitution variant inter-sample count correlations (all samples)",
      cut = as.factor(temp_dt[data.table::between(Nham_nt,1,2),Nham_nt]),
      size = 0.1,
      thresholds = list(
        "dotted" = log10(dimsum_meta[["fitnessMinInputCountAny"]] + 1), 
        "dashed" = log10(dimsum_meta[["fitnessMinInputCountAll"]] + 1)))
  }

  #Render report
  dimsum__render_report(dimsum_meta = dimsum_meta)

}

