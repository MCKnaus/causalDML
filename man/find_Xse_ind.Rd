% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plasso.R
\name{find_Xse_ind}
\alias{find_Xse_ind}
\title{Helper function finds the position for pre-specified SE rules}
\usage{
find_Xse_ind(CV, ind_min, oneSE, factor)
}
\arguments{
\item{CV}{Vector of cross-validated criterion}

\item{ind_min}{Index of cross-validated minimum}

\item{oneSE}{Standard error of cross-validated criterion at the minimum}

\item{factor}{Factor in which direction to go. Negative smaller model, positive larger model}
}
\value{
Index on the Lambda grid
}
\description{
Helper function finds the position for pre-specified SE rules
}
\keyword{internal}
