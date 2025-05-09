% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_metrics.R
\name{qc_metrics}
\alias{qc_metrics}
\title{Calculate standard QC metrics
This function calculates standard QC metrics.}
\usage{
qc_metrics(pipeline, mitochondrial_ens_ids)
}
\arguments{
\item{pipeline}{Kallisto or CellRanger SCE object}

\item{mitochondrial_ens_ids}{Mitochondrial genes (ENS_IDs)}
}
\value{
SCE object with quality metrics
}
\description{
Calculate standard QC metrics
This function calculates standard QC metrics.
}
