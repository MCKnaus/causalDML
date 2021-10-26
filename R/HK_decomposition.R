#' Heiler & Knaus decomposition of effect heterogeneity under non-homogeneous treatments
#'
#' This function builds on an \code{\link{APO_dml}} object to calculate the decomposition parameters
#' nATE(Z) = rATE(Z) + Delta(Z) described in Heiler and Knaus (2021).
#'
#' @param APO_dml \code{\link{APO_dml}} object containing the effective treatment specific doubly robust score
#' and nuisance parameters.
#' @param z Vector or matrix containing the low-dimensional heterogeneity variables to be considered. If \code{NULL},
#' the unconditional decomposition parameters are calculated.
#' @param intercept Adds an intercept to the variables in \code{z}. Set \code{FALSE} e.g. to get subgroup effects if
#' one-hot coded \code{z}.
#' @param spline If \code{TRUE}, B-spline regression based on \code{\link{crs}} is conducted.
#' @param subset Optional logical vector if decomposition should be run for subsample of original sample.
#' @param ... Pass option for \code{\link{crs}}.
#'
#' @return  \code{HK_decomposition} object with a list of nATE, rATE and Delta list containing each
#'  \item{score}{The parameter specific Neyman orthogonal score that was used as pseudo-outcome.}
#'  \item{list}{The \code{\link{lm}} object run to get the least squares coefficients (not with corrected standard errors).}
#'  \item{vcov}{The adjusted variance-covariance matrix derived in Heiler and Knaus (2021) for rATE and Delta.
#'  For nATE standard heteroskedasticity robust variance-covariance matrix.}
#'  \item{results}{A \code{\link{coeftest}} object with the coefficients and the appropriate inference results.}
#'
#' @references
#' \itemize{
#' \item Heiler, P., Knaus, M. C. (2021). Effect or Treatment Heterogeneity?
#' Policy Evaluation with Aggregated and Disaggregated Treatments
#' arXiv preprint arXiv:2110.01427. \url{https://arxiv.org/abs/2110.01427}
#' }
#'
#' @import lmtest
#'
#' @export
#'
HK_decomposition = function(APO_dml,z=NULL,intercept=TRUE,spline = FALSE,subset = NULL,...){

  # Extract relevant components of APO_dml object
  scores = APO_dml$gamma
  w = APO_dml$w_mat
  pscores = APO_dml$e_mat
  m_hat =  APO_dml$m_mat
  y =  APO_dml$y

  # Prepare heterogeneity variables
  if (is.null(z)) {
    z = matrix(1,dim(scores)[1],1)
    intercept = FALSE
  }
  else z = as.matrix(z)

  # Extract subset if specified
  if (!is.null(subset)) {
    scores = scores[subset,]
    w = w[subset,]
    pscores = pscores[subset,]
    m_hat = m_hat[subset,]
    y = y[subset]
    z = z[subset,]
  }

  # Define dimensions and additional variables
  n = dim(scores)[1]
  pi = colMeans(w)
  d = rowSums(w[,-1])

  # Initialize lists
  nATE = rATE = Delta = vector("list")

  # Get combined outcome nuisance for nATE
  m_hat_nATE = rowSums( m_hat[,-1] * (pscores[,-1]/rowSums(pscores[,-1])) )

  # Define decomposition scores
  nATE$score = m_hat_nATE + (y - m_hat_nATE) * d / rowSums(pscores[,-1]) - scores[,1]
  rATE$score = (scores[,-1] - scores[,1]) %*% (pi[-1]/sum(pi[-1]))
  Delta$score = nATE$score - rATE$score

  if (isTRUE(spline)) {
    # Run spline regressions
    fit_nATE = crs(nATE$score~z,...)
    fit_rATE = crs(rATE$score~z,...)
    fit_Delta = crs(Delta$score~z,...)

    # Keep the most flexible model as model matrix
    z = fit_nATE$model.lm$model
    if (ncol(fit_rATE$model.lm$model) > ncol(z)) z = fit_rATE$model.lm$model
    if (ncol(fit_Delta$model.lm$model) > ncol(z)) z = fit_Delta$model.lm$model
    z = z[,-1]
    intercept = TRUE
  }

  # Get Best linear Predictor
  if (isTRUE(intercept)){
    nATE$lm = lm(nATE$score ~ z)
    rATE$lm = lm(rATE$score ~ z)
    Delta$lm = lm(Delta$score ~ z)
  } else {
    nATE$lm = lm(nATE$score ~ 0 + z)
    rATE$lm = lm(rATE$score ~ 0 + z)
    Delta$lm = lm(Delta$score ~ 0 + z)
  }

  # Get Var-Cov matrix
  # Calculate adjustment term
  z = model.matrix(nATE$lm)
  k = ncol(z)
  Ebphi = matrix(0,k,length(pi)-1)
  ahat = matrix(0,n,k)
  for (k in seq(1,k)){
    Ebphi[k,] = colMeans( (scores[,-1] - scores[,1])*z[,k] )
    ahat[,k] = ( (w[,-1]*(1-pi[1]) + w[,1] %*% t(pi[-1])) %*% Ebphi[k,] ) / (1-pi[1])^2
  }

  ## Get vcov matrix
  # Get general components
  ZtZ = crossprod(z)
  Qinv = solve(ZtZ)

  # nATE
  ehat = nATE$lm$residuals
  be = z*ehat
  Sigma_nATE = crossprod(be)
  nATE$vcov = Qinv %*% Sigma_nATE %*% Qinv
  nATE$results = coeftest(nATE$lm,vcov. = nATE$vcov)

  # rATE
  ehat = rATE$lm$residuals
  be = z*ehat
  bea = be + sweep(ahat,2,colMeans(ahat))
  Sigma_rATE = crossprod(be)
  rATE$vcov = Qinv %*% Sigma_rATE %*% Qinv
  rATE$results = coeftest(rATE$lm,vcov. = rATE$vcov)

  # Delta
  ehat = Delta$lm$residuals
  be = z*ehat
  bea = be - sweep(ahat,2,colMeans(ahat))
  Sigma_Delta = crossprod(be)
  Delta$vcov = Qinv %*% Sigma_Delta %*% Qinv
  Delta$results = coeftest(Delta$lm,vcov. = Delta$vcov)

  output = list("nATE" = nATE, "rATE" = rATE, "Delta" = Delta)
  class(output) = "HK_decomposition"
  output
}


#' \code{summary} method for class \code{\link{HK_decomposition}}
#'
#' Prints the coefficients and inference of the three decomposition parameters.
#'
#' @param object \code{\link{HK_decomposition}}
#'
#' @import lmtest
#'
#' @export
#'
summary.HK_decomposition = function(object) {
  cat("nATE:")
  print(coeftest(object$nATE$lm,vcov. = object$nATE$vcov))
  cat("rATE:")
  print(coeftest(object$rATE$lm,vcov. = object$rATE$vcov))
  cat("Delta:")
  print(coeftest(object$Delta$lm,vcov. = object$Delta$vcov))
}


#' \code{predict} method for class \code{\link{HK_decomposition}}
#'
#' Provides the fitted values of the decomposition parameter of choice including adjusted standard error.
#'
#' @param object \code{\link{HK_decomposition}}
#' @param param chose from \code{c("nATE","rATE","Delta")}. Default \code{param="Delta"}
#'
#' @return n times 2 matrix with first column fitted value and second column standard error.
#'
#' @export
#'
predict.HK_decomposition = function(object,param="Delta") {
  z = model.matrix(object$nATE$lm)
  if (param == "nATE") {
    coef = object$nATE$lm$coefficients
    vcov = object$nATE$vcov
  } else if (param == "rATE") {
    coef = object$rATE$lm$coefficients
    vcov = object$rATE$vcov
  } else if (param == "Delta") {
    coef = object$Delta$lm$coefficients
    vcov = object$Delta$vcov
  } else stop("Choose first input: param from c(nATE, rATE, Delta) characters")

  pred = z %*% coef
  se = sqrt( rowSums((z %*% vcov) * z) )

  results = cbind(pred, se)
  return(results)
}


#' \code{plot} coefficients for class \code{\link{HK_decomposition}}
#'
#' @param object Object of class \code{\link{APO_dml}}.
#' @param sl Significance level for confidence intervals (default 0.05).
#'
#' @import ggplot2
#'
#' @export
#'
plot.HK_decomposition = function(object, sl=0.05) {
  group = Effect = cil = ciu = NULL

  qt = qt(1 - sl / 2, nobs(object$nATE$lm)-1)

  coef = c(object$nATE$lm$coefficients,
           object$rATE$lm$coefficients,
           object$Delta$lm$coefficients)
  se = c(sqrt(diag(object$nATE$vcov)),
         sqrt(diag(object$rATE$vcov)),
         sqrt(diag(object$Delta$vcov)))
  label_coef = colnames(model.matrix(object$nATE$lm))
  label_estimand = c("nATE","rATE","Delta")
  df = data.frame(Effect = coef,
             group = factor(rep(label_coef,3),levels = rev(label_coef)),
             cil = coef - qt * se,
             ciu = coef + qt * se,
             estimand = factor(rep(label_estimand,each=length(coef)/3),label_estimand))
  ggplot(df,aes(x=group,y=Effect,ymin=cil,ymax=ciu)) +
    geom_bar(stat="identity",color="darkgreen",fill="darkgreen",width = 0.5) +
    geom_errorbar(width = 0.3) +
    geom_hline(yintercept = 0) + theme_bw() + theme(axis.title.y = element_blank()) +
    facet_wrap(~estimand) + coord_flip()
}

