% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HK2_decomposition.R
\name{HK2_decomposition}
\alias{HK2_decomposition}
\title{Implementation of Heiler & Knaus (2025) decomposition.}
\usage{
HK2_decomposition(
  Y,
  A,
  G,
  T_mat,
  e_mat,
  m_mat,
  sampling_weights = NULL,
  cl = NULL
)
}
\arguments{
\item{Y}{Numeric vector containing the outcome variable.}

\item{A}{Aggregate treatment vector, e.g. binary. Program finds mapping to effective
treatment automatically if possible.}

\item{G}{Heterogeneity group vector. Provide as factor to control ordering,
otherwise program orders treatments in ascending order or alphabetically.}

\item{T_mat}{Logical matrix of effective treatment indicators (n x J). 
For example created by \code{\link{prep_w_mat}}.}

\item{e_mat}{n x J matrix with propensity scores.}

\item{m_mat}{n x J matrix with fitted outcome values.}

\item{sampling_weights}{Optional vector of sampling weights.}

\item{cl}{If not NULL, vector with cluster variables.}
}
\value{
Returns an \code{HK2_decomposition} object:
         \item{parameter}{14 x # of heterogeneity groups x # of treatment aggregates x 2 array
         storing point estimates and standard error of target and intermediate parameters. The ordering is
         c("CM","ACM","d0","s1","s2","s3","d1","d2","d3","Cov(etX,mut|Xg)","d4","SRCT2","d4'","d5").}
         \item{IFs}{14 x # of heterogeneity groups x # of treatment aggregates x n x 6 array
         storing the influence fcts and its components corresponding to each parameter for further use.}
         \item{mapping}{Logical matrix storing the mapping of aggregate and effective treatment.}
         \item{label}{List of labels for further useage.}
         \item{cl}{Cluster variable if specified.}
}
\description{
Implementation of Heiler & Knaus (2025) decomposition.
}
\references{
\itemize{
\item Heiler, P., Knaus, M.C. (2025). Heterogeneity analysis with heterogeneous treatments.
}
}
