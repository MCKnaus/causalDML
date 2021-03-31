#' Estimates non-paramteric CATEs using kernel regression as proposed by Fan et al. (2019)
#' and Zimmert & Lechner (2019).
#'
#' @param delta Vector of doubly robust ATE score. E.g obtained as one column of \code{ATE_dml$delta}
#' from \code{\link{ATE_dml}}.
#' @param z Heterogeneity variable(s) vector, matrix or data.frame.
#' @param bw_factor Factor by which cross-validated is multiplied. Default is undersmoothing with
#' factor 0.9 as recommended by Zimmert & Lechner (2019).
#'
#' @import np
#'
#' @return \code{kr_cate} object:
#' \item{model}{npqregression object (see\code{\link{npreg}}) of the kernel regression}
#' \item{fit}{Fitted values of the kernel regression}
#' \item{bw}{Cross-validated bandwidth (not scaled)}
#' \item{ate}{Average treatment effect}
#'
#' @references
#' \itemize{
#' \item Fan, Q., Hsu, Y.-C., Lieli, R. P., & Zhang, Y. (2019). Estimation of conditional average treatment effects with high-dimensional data.
#' arXiv preprint arXiv:1908.02399. \url{http://arxiv.org/abs/1908.02399}
#' \item Zimmert, M., & Lechner, M. (2019). Nonparametric estimation of causal heterogeneity under high-dimensional confounding.
#' arXiv preprint arXiv:1908.02399. \url{http://arxiv.org/abs/1908.08779}
#' }
#'
#' @export
#'
kr_cate = function(delta,z,bw_factor=0.9) {
  z = as.data.frame(z)
  delta = c(delta)
  bwobj = npregbw(ydat = delta, xdat = z, ckertype = 'gaussian', ckerorder = 2, regtype = 'lc', bwmethod = 'cv.ls')
  bws = bwobj$bw
  bw = bwobj$bw * bw_factor
  cate_model = npreg(tydat = delta, txdat = z, bws=bw, ckertype = 'gaussian', ckerorder = 2, regtype = 'lc')
  fit = fitted(cate_model)
  output = list("model"=cate_model,"fit"=fit,"bw"=bws,"ate"=mean(delta))
  class(output) = "kr_cate"
  output
}


#' \code{plot} method for class \code{\link{kr_cate}}
#'
#' @param kr_cate Object of class \code{\link{kr_cate}}.
#' @param z_label Optional label for the x-axis.
#' @param xrange Optional range of heterogeneity variable to be plotted on x-axis (e.g. \code{c(-1,1)}).
#' @param yrange Optional range of CATE to be plotted on y-axis (e.g. \code{c(-1,1)}).
#' @param sl Significance level for confidence intervals (default 0.05).
#'
#' @import ggplot2
#'
#' @export
#'
plot.kr_cate = function(kr_cate,
                        z_label=NULL,
                        xrange=NULL,
                        yrange=NULL,
                        sl=0.05) {
  z = cate = sehigh = selow = NULL

  qt = qt(1-sl / 2,length(kr_cate$fit)-1)
  df = data.frame(z = unlist(kr_cate$model$eval),cate = kr_cate$fit,
                  selow = kr_cate$fit - qt * se(kr_cate$model),
                  sehigh = kr_cate$fit + qt * se(kr_cate$model))
  if (!is.null(xrange)) df = df[df$z >= xrange[1] & df$z <= xrange[2],]

  g = ggplot(data = df,mapping = aes(x = z, y = cate)) +
    geom_line(size = .8) + geom_hline(yintercept=0) +
    geom_ribbon(aes(ymin = selow,max = sehigh),alpha=0.3) +
    theme_bw() + ylab("Conditional average treatment effect") +
    theme(text=element_text(size = 16, colour="black")) +
    theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") +
    geom_hline(yintercept=kr_cate$ate,linetype = 'dashed')

  if (!is.null(z_label)) g = g + xlab(z_label)
  if (!is.null(yrange)) g = g + coord_cartesian(ylim= yrange)
  return(g)
}


#' Estimates non-paramteric CATEs using spline regression as proposed by
#' Semenova and Chernozhukov (2020) using the \code{\link{crs}} package.
#'
#' @param delta Vector of doubly robust ATE score. E.g obtained as one column of \code{ATE_dml$delta}
#' from \code{\link{ATE_dml}}.
#' @param z Heterogeneity variable(s) vector, matrix or data.frame
#' @param ... Pass \code{\link{crs}} options, otherwise default settings are used
#'
#' @import crs
#'
#' @return \code{spline_cate} object:
#' \item{model}{\code{\link{crs}} object of the spline regression}
#' \item{fit}{Fitted values of the spline regression}
#' \item{se_fit}{Robust standard errors of fitted values}
#' \item{ate}{Average teratment effect}
#' \item{z}{Original heterogeneity variable(s) data.frame}
#'
#' @references
#' \itemize{
#' \item Semenova, V., & Chernozhukov, V. (2020). Debiased machine learning of conditional
#' average treatment effects and other causal functions. The Econometrics Journal,
#' utaa027, doi: https://doi.org/10.1093/ectj/utaa027
#' }
#'
#' @export
#'
spline_cate = function(delta,z, ...) {
  # Run crs
  z = as.matrix(z,nrow = length(delta))
  cate_model = crs(delta ~ z, ...)
  fit = cate_model$fitted.values
  # Get robust standard errors for all fitted values
  bz = as.matrix(cate_model$model.lm$model[,-1])
  if (ncol(bz) < cate_model$model.lm$rank) bz = cbind( rep(1,nrow(bz)), bz )
  bztbz = crossprod(bz)
  bztbz_inv = solve(bztbz)
  AV = bztbz_inv %*% crossprod(bz * cate_model$model.lm$residuals) %*% bztbz_inv
  # se_fit = sqrt( diag(bz %*% AV %*% t(bz)) ) https://stackoverflow.com/questions/21708489/compute-only-diagonals-of-matrix-multiplication-in-r
  se_fit = sqrt( rowSums((bz %*% AV) * bz) )

  output = list("model"=cate_model,"fit"=fit,"se_fit"=se_fit,"ate"=mean(delta),"z"=z)
  class(output) = "spline_cate"
  output
}


#' \code{plot} method for class \code{\link{spline_cate}}
#'
#' @param spline_cate Object of class \code{\link{spline_cate}}.
#' @param z_label Optional label for the x-axis.
#' @param xrange Optional range of heterogeneity variable to be plotted on x-axis (e.g. \code{c(-1,1)}).
#' @param yrange Optional range of CATE to be plotted on y-axis (e.g. \code{c(-1,1)}).
#' @param sl Significance level for confidence intervals (default 0.05).
#'
#' @import ggplot2
#'
#' @export
#'
plot.spline_cate = function(spline_cate,
                            z_label=NULL,
                            xrange=NULL,
                            yrange=NULL,
                            sl=0.05) {
  if (ncol(spline_cate$z) > 1) stop("Currently only supports one heterogeneity variable.")

  qt = qt(1-sl / 2,length(spline_cate$fit)-1)
  df = data.frame(z = spline_cate$z,cate = spline_cate$fit,
                  selow = spline_cate$fit - qt * spline_cate$se_fit,
                  sehigh = spline_cate$fit + qt * spline_cate$se_fit)
  if (!is.null(xrange)) df = df[df$z >= xrange[1] & df$z <= xrange[2],]

  g = ggplot(data = df,mapping = aes(x = z, y = cate)) +
    geom_line(size = .8) + geom_hline(yintercept=0) +
    geom_ribbon(aes(ymin = selow,max = sehigh),alpha=0.3) +
    theme_bw() + ylab("Conditional average treatment effect") +
    theme(text=element_text(size = 16, colour="black")) +
    theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") +
    geom_hline(yintercept=spline_cate$ate,linetype = 'dashed')

  if (!is.null(z_label)) g = g + xlab(z_label)
  if (!is.null(yrange)) g = g + coord_cartesian(ylim= yrange)
  return(g)
}


