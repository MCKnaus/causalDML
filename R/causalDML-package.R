#' causalDML: Causal Double Machine Learning
#'
#' An implementation of Double Machine Learning based method
#' as they are discussed by Chernozhukov et al. (2018) and reviewed in Knaus (2020).
#' It is tailored for program evaluation with potentially multiple treatments estimates
#' average potential outcomes and average treatment effects.
#'
#' @docType package
#'
#' @author Michael C. Knaus
#'
#' @examples
#' # Generate data
#' n = 50
#' p = 10
#' X = matrix(rnorm(n * p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#'
#' # Run causal DML
#' cDML = causalDML(Y, W, X)
#'
#' # Show average potential outcome estimates
#' summary(cDML$APO)
#'
#' # Show average tretment effect estimates
#' summary(cDML$ATE)
#'
#' @references
#' \itemize{
#' \item Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., &Robins, J. (2018).
#' Double/Debiased machine learning for treatment and structuralparameters.The Econometrics Journal,21(1), C1-C68
#' \item Knaus, M. C. (2020). Double machine learning based program evaluation under unconfoundedness.
#'   arXiv preprint arXiv:2003.03191.\url{http://arxiv.org/abs/2003.03191}
#' }
#'
#' @name causalDML
#'
#' @useDynLib causalDML, .registration = TRUE
NULL
