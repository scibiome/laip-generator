% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LAIP_Generator.R
\name{print_scatterplot}
\alias{print_scatterplot}
\title{print_scatterplot
print_scatterplot is a function for debugging purpose}
\usage{
print_scatterplot(
  data,
  max_rows = 10000,
  x = "SSC-H",
  y,
  log_transform_y = TRUE,
  laipColumnName = NULL
)
}
\arguments{
\item{data}{data.frame or matrix with events to be plotted}

\item{max_rows}{numeric (default = 10000), specifying the number of events to be plotted}

\item{x}{string, specifying the parameter to be plotted on x-axis}

\item{y}{string, specifying the parameter to be plotted on y-axis}

\item{log_transform_y}{logical (default: TRUE); if TRUE the values to be plotted on y-axis are log_transformed}

\item{laipColumnName}{string, specifying the parameter marking simulated LAIP-positive events}
}
\value{

}
\description{
print_scatterplot
print_scatterplot is a function for debugging purpose
}
