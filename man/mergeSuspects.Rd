% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/masses.R
\name{mergeSuspects}
\alias{mergeSuspects}
\title{Merge suspect lists}
\usage{
mergeSuspects(massSuspects, formulaSuspects)
}
\arguments{
\item{massSuspects}{data frame with `name`,`mass` columns}

\item{formulaSuspects}{data frame with `name`,`mass` columns}
}
\description{
Merges suspect lists generated from formula calculations and from mass calculations.
If formula results are present, they take priority over mass results. Note that formula results must first be
converted to masses using [convertFormulas].
}
\author{
stravsmi
}
