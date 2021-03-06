% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoFRA.R
\name{filterData}
\alias{filterData}
\title{This function filter a data frame on column named "pro"}
\usage{
filterData(dfPro, filter)
}
\arguments{
\item{dfPro}{data frame with iBAQ values}

\item{filter}{character vector with items to remove}
}
\description{
This function filter a data frame on column named "pro"
}
\examples{
library(CoFRA)
data(iBAQ)
iBAQ2=filterData(iBAQ,">CON") # filter headers starting with >CON
}
\keyword{filter}

