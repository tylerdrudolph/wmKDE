% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcHR.R
\name{calcHR}
\alias{calcHR}
\title{calcHR}
\usage{
calcHR(est, p = c(0.5, 0.95), silent = TRUE)
}
\arguments{
\item{est}{UD object, list containing 3 elements: x (x1) and y (x2) coordinates and corresponding density estimates (fhat)}

\item{p}{vector of probabilities at which to derive home range size}

\item{silent}{logical indicating whether or not to print messages (default = TRUE)}
}
\value{
probability density bin values associated with user-defined isopleths
}
\description{
Function to calculate home range size and isopleth values corresponding to UD probability contours. Adapted from code shared by John Fieberg.
}
