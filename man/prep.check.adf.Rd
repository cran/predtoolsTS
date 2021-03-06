% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep.R
\name{prep.check.adf}
\alias{prep.check.adf}
\title{Augmented Dickey-Fuller test}
\usage{
prep.check.adf(tserie)
}
\arguments{
\item{tserie}{a \code{ts} or a \code{prep} object}
}
\description{
Performs ADF test just as another tool to check stationarity.
}
\details{
Shows the results of an ADF test. A p-value<0.05 suggests the data is stationary.
}
\examples{
prep.check.adf(AirPassengers)
prep.check.adf(prep(AirPassengers))
}
