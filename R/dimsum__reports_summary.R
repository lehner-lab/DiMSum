
#' dimsum__reports_summary
#'
#' Generate a summary report in html format.
#'
#' @param dimsum_meta an experiment metadata object (required)
#'
#' @return a summary report character string (html format)
#' @export
dimsum__reports_summary <- function(dimsum_meta){
  #HTML results summary
  #Piece1
  reports_piece1 <- "
   <html>\
   <head>\
    <title>PROJECT_NAME</title>\
   </head>\
   <body>\
   \t<h1>DiMSum DIMSUM_VERSION Pipeline Report (PROJECT_NAME)</h1>\
   \t\
   \t<div>\
   \t\t<a href=\"#FASTQC\"> 1. Raw FASTQ file quality control (FASTQC) </a><br>\
   \t\t<a href=\"#CUTADAPT\"> 2. Trimming of constant regions (cutadapt) </a><br>\
   \t\t<a href=\"#USEARCH1\"> 3. Alignment of overlapping paired-end reads (USEARCH) </a><br>\
   \t\t<a href=\"#USEARCH2\"> 4. Aligned length (USEARCH) </a><br>\
   \t\t<a href=\"#MERGE1\"> 5. Nucleotide mutation statistics </a><br>\
   \t\t<a href=\"#MERGE2\"> 6. Amino acid mutation statistics </a><br>\
   \t\t<a href=\"#COUNTDISTR\"> 7. Variant count distributions </a><br>\
   \t\t<a href=\"#COUNTSCATTERMATRIX\"> 8. Inter-sample variant count diagnostic plot </a><br>\
   \t\t<a href=\"#INPUTTHRESHOLD\"> 9. Input read threshold for full fitness range </a><br>\
   \t\t<a href=\"#REPLICATEDEVIATION\"> 10. Affect of fitness normalisation on replicate deviations </a><br>\
   \t\t<a href=\"#FITNESSDISTR\"> 11. Fitness distributions </a><br>\
   \t\t<a href=\"#FITNESSSCATTERMATRIX\"> 12. Inter-sample variant fitness diagnostic plot </a><br>\
   \t\t<a href=\"#FITNESSERROR\"> 13. Fitness error model </a><br>\
   \t</div>\
  \
   \t<div id=\"FASTQC\">\
   \t\t<h2>1. Raw FASTQ file quality control (FASTQC)</h2>\
   \t\t<p><img src=\"reports/dimsum_stage_fastqc_report_pair1_fastqc.png\" width=\"1200\"></p>\
  "
  #Piece1 (nucleotide sequence)
  reports_piece1_nt <- "
   <html>\
   <head>\
    <title>PROJECT_NAME</title>\
   </head>\
   <body>\
   \t<h1>DiMSum DIMSUM_VERSION Pipeline Report (PROJECT_NAME)</h1>\
   \t\
   \t<div>\
   \t\t<a href=\"#FASTQC\"> 1. Raw FASTQ file quality control (FASTQC) </a><br>\
   \t\t<a href=\"#CUTADAPT\"> 2. Trimming of constant regions (cutadapt) </a><br>\
   \t\t<a href=\"#USEARCH1\"> 3. Alignment of overlapping paired-end reads (USEARCH) </a><br>\
   \t\t<a href=\"#USEARCH2\"> 4. Aligned length (USEARCH) </a><br>\
   \t\t<a href=\"#MERGE1\"> 5. Nucleotide mutation statistics </a><br>\
   \t\t<a href=\"#COUNTDISTR\"> 6. Variant count distributions </a><br>\
   \t\t<a href=\"#SCATTERMATRIX\"> 7. Inter-sample variant count diagnostic plot </a><br>\
   \t\t<a href=\"#INPUTTHRESHOLD\"> 8. Input read threshold for full fitness range </a><br>\
   \t\t<a href=\"#REPLICATEDEVIATION\"> 9. Affect of fitness normalisation on replicate deviations </a><br>\
   \t\t<a href=\"#FITNESSDISTR\"> 10. Fitness distributions </a><br>\
   \t\t<a href=\"#FITNESSSCATTERMATRIX\"> 11. Inter-sample variant fitness diagnostic plot </a><br>\
   \t\t<a href=\"#FITNESSERROR\"> 12. Fitness error model </a><br>\
   \t</div>\
  \
   \t<div id=\"FASTQC\">\
   \t\t<h2>1. Raw FASTQ file quality control (FASTQC)</h2>\
   \t\t<p><img src=\"reports/dimsum_stage_fastqc_report_pair1_fastqc.png\" width=\"1200\"></p>\
  "
  #Paired-end read piece1
  reports_PE_piece1 <- "
   \t\t<p><img src=\"reports/dimsum_stage_fastqc_report_pair2_fastqc.png\" width=\"1200\"></p>\
  "
  #Piece2
  reports_piece2 <- "
   \t</div>\
  \
    \t<div id=\"CUTADAPT\">\
    \t\t<h2>2. Trimming of constant regions (cutadapt)</h2>\
  \t \t<p><img src=\"reports/dimsum_stage_cutadapt_report_pair1.png\" width=\"1200\"></p>\
  "
  #Paired-end read piece2
  reports_PE_piece2 <- "
  \t \t<p><img src=\"reports/dimsum_stage_cutadapt_report_pair2.png\" width=\"1200\"></p>\
  "
  #Piece3
  reports_piece3 <- "
   \t</div>\
   \t\
    \t<div id=\"USEARCH1\">\
    \t\t<h2>3. Alignment of overlapping paired-end reads (USEARCH)</h2>\
  \t \t<p><img src=\"reports/dimsum_stage_usearch_report_paircounts.png\" width=\"1200\"></p>\
   \t</div>\
  \
    \t<div id=\"USEARCH2\">\
    \t\t<h2>4. Aligned length (USEARCH)</h2>\
  \t \t<p><img src=\"reports/dimsum_stage_usearch_report_mergedlength.png\" width=\"1200\"></p>\
   \t</div>\
  \
      <div id=\"MERGE1\">\
        <h2>5. Nucleotide mutation statistics</h2>\
      <p><img src=\"reports/dimsum_stage_merge_report_nucmutationpercentages.png\" width=\"1200\"></p>\
      <p><img src=\"reports/dimsum_stage_merge_report_nucmutationcounts.png\" width=\"1200\"></p>\
    </div>\
  \
    \t<div id=\"MERGE2\">\
    \t\t<h2>6. Amino acid mutation statistics</h2>\
  \t \t<p><img src=\"reports/dimsum_stage_merge_report_aamutationpercentages.png\" width=\"1200\"></p>\
  \t \t<p><img src=\"reports/dimsum_stage_merge_report_aamutationcounts.png\" width=\"1200\"></p>\
   \t</div>\
  \
      <div id=\"COUNTDISTR\">\
        <h2>7. Variant count distributions</h2>\
      <p><img src=\"reports/dimsum_stage_diagnostics_report_count_hist_input.png\" width=\"1200\"></p>\
      <p><img src=\"reports/dimsum_stage_diagnostics_report_count_hist_output.png\" width=\"1200\"></p>\
    </div>\
  \
      <div id=\"COUNTSCATTERMATRIX\">\
        <h2>8. Inter-sample variant count diagnostic plot</h2>\
      <p><img src=\"reports/dimsum_stage_diagnostics_report_scatterplotmatrix_all.png\" width=\"1200\"></p>\
    </div>\
  \
      <div id=\"INPUTTHRESHOLD\">\
        <h2>9. Input read threshold for full fitness range</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.png\" width=\"600\"></p>\
    </div>\
  \
      <div id=\"REPLICATEDEVIATION\">\
        <h2>10. Affect of fitness normalisation on replicate deviations</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicate_deviation_scatter.png\" width=\"600\"></p>\
    </div>\
  \
      <div id=\"FITNESSDISTR\">\
        <h2>11. Fitness distributions (before and after inter-replicate normalisation)</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.png\" width=\"600\"><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.png\" width=\"600\"></p>\
    </div>\
  \
      <div id=\"FITNESSSCATTERMATRIX\">\
        <h2>12. Inter-sample variant fitness diagnostic plot (before and after inter-replicate normalisation)</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.png\" width=\"1000\"></p>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.png\" width=\"1000\"></p>\
    </div>\
  \
      <div id=\"FITNESSERROR\">\
        <h2>13. Fitness error model</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_repspec.png\" width=\"800\"></p>\
    </div>\
  \
   </body>\
  </html>"
  #Piece3 (nucleotide sequence)
  reports_piece3_nt <- "
   \t</div>\
   \t\
    \t<div id=\"USEARCH1\">\
    \t\t<h2>3. Alignment of overlapping paired-end reads (USEARCH)</h2>\
  \t \t<p><img src=\"reports/dimsum_stage_usearch_report_paircounts.png\" width=\"1200\"></p>\
   \t</div>\
  \
    \t<div id=\"USEARCH2\">\
    \t\t<h2>4. Aligned length (USEARCH)</h2>\
  \t \t<p><img src=\"reports/dimsum_stage_usearch_report_mergedlength.png\" width=\"1200\"></p>\
   \t</div>\
  \
      <div id=\"MERGE1\">\
        <h2>5. Nucleotide mutation statistics</h2>\
      <p><img src=\"reports/dimsum_stage_merge_report_nucmutationpercentages.png\" width=\"1200\"></p>\
      <p><img src=\"reports/dimsum_stage_merge_report_nucmutationcounts.png\" width=\"1200\"></p>\
    </div>\
  \
      <div id=\"COUNTDISTR\">\
        <h2>6. Variant count distributions</h2>\
      <p><img src=\"reports/dimsum_stage_diagnostics_report_count_hist_input.png\" width=\"1200\"></p>\
      <p><img src=\"reports/dimsum_stage_diagnostics_report_count_hist_output.png\" width=\"1200\"></p>\
    </div>\
  \
      <div id=\"COUNTSCATTERMATRIX\">\
        <h2>7. Inter-sample variant count diagnostic plot</h2>\
      <p><img src=\"reports/dimsum_stage_diagnostics_report_scatterplotmatrix_all.png\" width=\"1200\"></p>\
    </div>\
  \
      <div id=\"INPUTTHRESHOLD\">\
        <h2>8. Input read threshold for full fitness range</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_inputcounts.png\" width=\"600\"></p>\
    </div>\
  \
      <div id=\"REPLICATEDEVIATION\">\
        <h2>9. Affect of fitness normalisation on replicate deviations</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicate_deviation_scatter.png\" width=\"600\"></p>\
    </div>\
  \
      <div id=\"FITNESSDISTR\">\
        <h2>10. Fitness distributions (before and after inter-replicate normalisation)</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density.png\" width=\"600\"><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicates_density_norm.png\" width=\"600\"></p>\
    </div>\
  \
      <div id=\"FITNESSSCATTERMATRIX\">\
        <h2>11. Inter-sample variant fitness diagnostic plot (before and after inter-replicate normalisation)</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter.png\" width=\"1000\"></p>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_fitness_replicates_scatter_norm.png\" width=\"1000\"></p>\
    </div>\
  \
      <div id=\"FITNESSERROR\">\
        <h2>12. Fitness error model</h2>\
      <p><img src=\"reports/dimsum_stage_fitness_report_1_errormodel_repspec.png\" width=\"800\"></p>\
    </div>\
  \
   </body>\
  </html>"
  #Assemble HTML
  reports_summary_template <- ""
  temp_reports_piece1 <- ""
  temp_reports_piece3 <- ""
  if(dimsum_meta[["sequenceType"]]=="coding"){
    temp_reports_piece1 <- reports_piece1
    temp_reports_piece3 <- reports_piece3
  }else{
    temp_reports_piece1 <- reports_piece1_nt
    temp_reports_piece3 <- reports_piece3_nt
  }
  if(dimsum_meta[["paired"]]){
    reports_summary_template <- paste0(temp_reports_piece1, reports_PE_piece1, reports_piece2, reports_PE_piece2, temp_reports_piece3)
  }else{
    reports_summary_template <- paste0(temp_reports_piece1, reports_piece2, temp_reports_piece3)
  }
  reports_summary_template <- gsub("PROJECT_NAME", dimsum_meta[["projectName"]], reports_summary_template)
  reports_summary_template <- gsub("DIMSUM_VERSION", as.character(packageVersion("DiMSum")), reports_summary_template)
  return(reports_summary_template)
}
