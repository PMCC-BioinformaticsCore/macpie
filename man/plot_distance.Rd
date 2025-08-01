% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_distance.R
\name{plot_distance}
\alias{plot_distance}
\title{Create a distance heatmap}
\usage{
plot_distance(data = NULL, group_by = NULL, treatment = NULL)
}
\arguments{
\item{data}{A Seurat object}

\item{group_by}{A metadata column name to group data}

\item{treatment}{to specify one treatment group in the group_by parameter}
}
\value{
it returns a pheatmap object
}
\description{
Plot a heatmap to show poisson distance on count matrix of pre-processed data.
The poisson distance matrix is calculated on subset of wells that specify by two
other parameters. By specifying group_by to a column in the metadata and a specific
treatment of interest, poisson distance is then calculated for corresponding wells.
}
\examples{
data(mini_mac)
p <- plot_distance(mini_mac, group_by = "combined_id", treatment = "DMSO_0")
}
