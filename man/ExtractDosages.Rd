% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ExtractDosages}
\alias{ExtractDosages}
\title{Function to extract SNPs from a binary dosage file}
\usage{
ExtractDosages(bdosageFilename, mapFilename, numSub, numSNPs)
}
\arguments{
\item{bdosageFilename}{Name of binary dosage file}

\item{mapFilename}{Name of map file associated with dosage file}

\item{numSub}{Number of subjects with data in dosage file}

\item{numSNPs}{Number of SNPs to read in}
}
\value{
List with a vector of dosages and a matrix of probabilities
and a list of the input values
}
\description{
Function to extract SNPs from a binary dosage file
}
