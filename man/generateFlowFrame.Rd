% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LAIP_Generator.R
\name{generateFlowFrame}
\alias{generateFlowFrame}
\title{generateFlowFrame}
\usage{
generateFlowFrame(data, parameters = NULL, description = NULL)
}
\arguments{
\item{data}{= data.frame containing the measured intensities. Rows correspond to cells, columns to the different measurement channels.}

\item{parameters}{= data.frame containing information about each column of the flowFrame. columns of parameters: name, description.}

\item{description}{= a list containing the meta data included in the FCS file.}
}
\value{
flowFrame-object
}
\description{
Generates a flowFrame based on the given data.frame, parameters and description.
Parameters and description are optional. Without given description the first
7 columns must be non-color-parameter (time, SSC-H, etc.) data
}
\examples{
generateFlowFrame(bm_data)
}
