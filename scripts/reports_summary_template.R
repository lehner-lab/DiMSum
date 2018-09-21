
#HTML results summary
reports_summary_template <- "<html>\
 <head>\
  <title>PROJECT_NAME</title>\
 </head>\
 <body>\
 \t<h1>DiMSum v0.1 Pipeline Report (PROJECT_NAME)</h1>\
 \t\
 \t<div>\
 \t\t<a href=\"#FASTQC\"> 1. Raw FASTQ file quality control (FASTQC) </a><br>\
 \t\t<a href=\"#CUTADAPT\"> 2. Trimming of constant regions (cutadapt) </a><br>\
 \t\t<a href=\"#USEARCH1\"> 3. Alignment of overlapping paired-end reads (USEARCH) </a><br>\
 \t\t<a href=\"#USEARCH2\"> 4. Alignment length (USEARCH) </a><br>\
 \t\t<a href=\"#MERGE1\"> 5. Nucleotide mutation statistics </a><br>\
 \t\t<a href=\"#MERGE2\"> 6. Amino acid mutation statistics </a><br>\
 \t\t<a href=\"#SUMMARY\"> 7. Variant filtering summary statistics </a><br>\
 \t\t<a href=\"#DIAGNOSTIC\"> 8. Variant count diagnostics </a><br>\
 \t</div>\
\
 \t<div id=\"FASTQC\">\
 \t\t<h2>1. Raw FASTQ file quality control (FASTQC)</h2>\
 \t\t<p><img src=\"reports/dimsum_stage_fastqc_report_pair1_fastqc.png\" width=\"1200\"></p>\
 \t\t<p><img src=\"reports/dimsum_stage_fastqc_report_pair2_fastqc.png\" width=\"1200\"></p>\
 \t</div>\
\
  \t<div id=\"CUTADAPT\">\
  \t\t<h2>2. Trimming of constant regions (cutadapt)</h2>\
\t \t<p><img src=\"reports/dimsum_stage_cutadapt_report_pair1.png\" width=\"1200\"></p>\
\t \t<p><img src=\"reports/dimsum_stage_cutadapt_report_pair2.png\" width=\"1200\"></p>\
 \t</div>\
 \t\
  \t<div id=\"USEARCH1\">\
  \t\t<h2>3. Alignment of overlapping paired-end reads (USEARCH)</h2>\
\t \t<p><img src=\"reports/dimsum_stage_usearch_report_paircounts.png\" width=\"1200\"></p>\
 \t</div>\
\
  \t<div id=\"USEARCH2\">\
  \t\t<h2>4. Alignment length (USEARCH)</h2>\
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
    <div id=\"SUMMARY\">\
      <h2>7. Variant filtering statistics</h2>\
    <p><img src=\"reports/dimsum_stage_merge_report_variantpercentages.png\" width=\"1200\"></p>\
    <p><img src=\"reports/dimsum_stage_merge_report_variantcounts.png\" width=\"1200\"></p>\
  </div>\
\
    <div id=\"DIAGNOSTIC\">\
      <h2>8. Variant count diagnostics</h2>\
    <p><img src=\"reports/dimsum_stage_diagnostics_report_count_hist_input.png\" width=\"1200\"></p>\
    <p><img src=\"reports/dimsum_stage_diagnostics_report_count_hist_output.png\" width=\"1200\"></p>\
    <p><img src=\"reports/dimsum_stage_diagnostics_report_scatterplotmatrix_input.png\" width=\"1200\"></p>\
    <p><img src=\"reports/dimsum_stage_diagnostics_report_scatterplotmatrix_output.png\" width=\"1200\"></p>\
  </div>\
\
 </body>\
</html>"
