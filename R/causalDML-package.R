#' causalDML: Causal Double Machine Learning
#'
#' An implementation of Double Machine Learning based method
#' as discussed by Chernozhukov et al. (2018) and reviewed in Knaus (2020).
#' Tailored for program evaluation with potentially multiple treatments, it estimates
#' average potential outcomes and average treatment effects.
#'
#' @keywords internal
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
#' # Show average treatment effect estimates
#' summary(cDML$ATE)
#'
#' @references
#' \itemize{
#' \item Chernozhukov, V., et al. (2018). Double/Debiased machine learning for treatment and structural parameters.
#'   *The Econometrics Journal*, 21(1), C1–C68.
#' \item Knaus, M. C. (2022). Double machine learning based program evaluation under unconfoundedness. 
#' *The Econometrics Journal*, 25(3). \url{https://doi.org/10.1093/ectj/utac015}
#' }
#'
#' @useDynLib causalDML, .registration = TRUE
"_PACKAGE"
