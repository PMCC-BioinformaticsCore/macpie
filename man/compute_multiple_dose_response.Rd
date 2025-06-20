% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_multiple_dose_response.R
\name{compute_multiple_dose_response}
\alias{compute_multiple_dose_response}
\title{Model Gene Dose-Response Curves Across Treatments}
\usage{
compute_multiple_dose_response(
  data,
  genes = NULL,
  normalisation = "limma_voom",
  control_value = "DMSO",
  batch = 1,
  k = 2,
  num_cores = 1
)
}
\arguments{
\item{data}{A Seurat or TidySeurat object containing expression data and metadata.}

\item{genes}{A character vector of gene names to model. If NULL, all significant DE genes across comparisons are used.}

\item{normalisation}{A character string indicating the normalization method. One of: "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
"edgeR", "RUVg", "RUVs", "RUVr", "limma_voom", "zinb". Default is "limma_voom".}

\item{control_value}{A string indicating the control condition in "Treatment_1". Default is "DMSO".}

\item{batch}{Batch variable to use for normalization if applicable. Default is 1.}

\item{k}{Number of unwanted factors for RUV normalization. Default is 2.}

\item{num_cores}{Number of CPU cores to use in parallel model fitting. Default is 1.}
}
\value{
A data frame of EC50 values per gene and treatment.
}
\description{
This function fits dose-response models for a set of genes across different treatments
using the \code{drc} package. It returns EC50 values per gene per treatment.
}
\examples{
\dontrun{
rds_file <- system.file("extdata/PMMSq033/PMMSq033.rds", package = "macpie")
mac <- readRDS(rds_file)
res <- compute_multiple_dose_response(
  data = mac,
  genes = c("PTPRA", "MYC"),
  normalisation = "limma_voom",
  treatment_value = "Camptothecin",
  num_cores = 2
)
head(res)
}
}
