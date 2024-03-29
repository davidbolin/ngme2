% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{ngme_ts_make_A}
\alias{ngme_ts_make_A}
\title{Make observation matrix for time series}
\usage{
ngme_ts_make_A(loc, replicate = NULL, range = c(min(loc), max(loc)))
}
\arguments{
\item{loc}{integers (after sorting, no gaps > 1)}

\item{replicate}{indicating replicate measure at same location}

\item{range}{range for the mesh
by default range=(min(loc), max(loc))}
}
\value{
A matrix (length(loc) * length(unique(loc)))
}
\description{
Make observation matrix for time series
}
\examples{
ngme_ts_make_A(c(1, 2, 2), replicate = c(1, 1, 2))
ngme_ts_make_A(c(1, 2, 2), range = c(1, 5))
}
