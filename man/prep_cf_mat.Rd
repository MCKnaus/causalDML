% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_effects.R
\name{prep_cf_mat}
\alias{prep_cf_mat}
\title{Creates matrix of binary cross-fitting fold indicators (n x # cross-folds)}
\usage{
prep_cf_mat(n, cf, w_mat = NULL, cl = NULL)
}
\arguments{
\item{n}{Number of observations.}

\item{cf}{Number of cross-fitting folds.}

\item{w_mat}{Optional logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
If specified, cross-fitting folds will preserve the treatment rations from full sample.}

\item{cl}{Optional vector of cluster variable if cross-fitting should account for clusters.}
}
\value{
Logical matrix of cross-fitting folds (n x # folds).
}
\description{
Creates matrix of binary cross-fitting fold indicators (n x # cross-folds)
}
