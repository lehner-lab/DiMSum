% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimsum__swap_reads.R
\name{dimsum__swap_reads}
\alias{dimsum__swap_reads}
\title{dimsum__swap_reads}
\usage{
dimsum__swap_reads(dimsum_meta, cutadapt_outpath)
}
\arguments{
\item{dimsum_meta}{an experiment metadata object (required)}

\item{cutadapt_outpath}{cutadapt output path (required)}
}
\value{
An updated dimsum_meta object
}
\description{
Swap reads in fastq files according to adapter presence.
}
