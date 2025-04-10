% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HK2_decomposition.R
\name{plot.HK2_decomposition}
\alias{plot.HK2_decomposition}
\title{\code{plot} method for class \code{\link{HK2_decomposition}}}
\usage{
\method{plot}{HK2_decomposition}(
  HK2_decomposition,
  decomposition = "adim",
  t_aggregate = c(3, 2),
  x_aggregate = c(3, 2),
  levels = F,
  pe_digits = 3
)
}
\arguments{
\item{HK2_decomposition}{Object of class \code{\link{HK2_decomposition}}.}

\item{decomposition}{Either "adim" or "dim" for the respective decomposition. Default "adim".}

\item{t_aggregate}{Integer scalar of vector of length two specifying treatment aggregate to be collected. 
1 is reserved for collection of all treaments. 
If single integer, level is extracted. 
Vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.}

\item{x_aggregate}{Integer scalar of vector of length two specifying covariate aggregate to be collected. 
If single integer, level is extracted. 
1 is reserved for the unconditional version. 
vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.}

\item{levels}{If TRUE, prints levels next to difference. Only applicable with length(x_aggregate) == 2.}

\item{pe_digits}{Controls number of digits printed for point estimate.}
}
\description{
\code{plot} method for class \code{\link{HK2_decomposition}}
}
