% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_multi_screen_profile.R
\name{plot_multi_screen_profile}
\alias{plot_multi_screen_profile}
\title{Plot multi-screen profile from fgsea results}
\usage{
plot_multi_screen_profile(
  data,
  color_by = NULL,
  size_by = "logPadj",
  label_angle = 90
)
}
\arguments{
\item{data}{A tidyseurat object with \code{screen_profile} in the @tools slot.}

\item{color_by}{A string specifying which column in screen_profile to use for color (default: automatically chosen).}

\item{size_by}{A string specifying which column in screen_profile to use for point size (default: log10(padj)).}

\item{label_angle}{Angle of x-axis labels (default: 90).}
}
\value{
A ggplot object
}
\description{
Plot multi-screen profile from fgsea results
}
\examples{
data(mini_mac)
mini_mac <- compute_multi_screen_profile(mini_mac, target = "Staurosporine_10")
plot_multi_screen_profile(mini_mac, color_by = "Sample_type")
}
