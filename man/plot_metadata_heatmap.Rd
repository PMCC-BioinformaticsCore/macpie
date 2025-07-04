% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_metadata_heatmap.R
\name{plot_metadata_heatmap}
\alias{plot_metadata_heatmap}
\title{Generate Heatmaps of Metadata Function}
\usage{
plot_metadata_heatmap(
  metadata = NULL,
  metadata_file = NULL,
  legend = TRUE,
  output_file = NULL,
  plate = NULL
)
}
\arguments{
\item{metadata}{Metadata object of the class data frame, matrix or tibble.}

\item{metadata_file}{Path to the metadata CSV file path.}

\item{legend}{A character value ("show" or "none") to control whether to display legends.}

\item{output_file}{A file path and name to save the heatmaps as a graph (png, pdf, or jpg).}

\item{plate}{Plate to be plotted. Default is first plate.}
}
\value{
Displays the plot as a ggplot object or saves it as a file (jpg, png, or pdf).
}
\description{
This function generates heatmaps from metadata, either from an existing
object or from a CSV file.
}
\examples{
#Example
metadata_file_path <- system.file("extdata", "PMMSq033_metadata.csv", package = "macpie")
metadata<-read_metadata(metadata_file_path)
plot_metadata_heatmap(metadata=metadata)
}
