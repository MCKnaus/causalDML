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
#' \item{ate}{Average teratment effect}
#'
#' @references
#' \itemize{
#' \item Fan, Q., Hsu, Y.-C., Lieli, R. P., & Zhang, Y. (2019). Estimation of conditional averagetreatment effects with high-dimensional data.
#' arXiv preprint arXiv:1908.02399. \url{http://arxiv.org/abs/1908.02399}
#' \item Zimmert, M., & Lechner, M. (2019). Nonparametric estimation of causal heterogeneityunder high-dimensional confounding.
#' arXiv preprint arXiv:1908.02399. \url{http://arxiv.org/abs/1908.08779}
#' }
#'
#' @export
#'
kr_cate = function(delta,z,bw_factor=0.9) {
  z = as.data.frame(z)
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
    theme(text=element_text(family="serif",size = 16, colour="black")) +
    theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") +
    geom_hline(yintercept=kr_cate$ate,linetype = 'dashed')

  if (!is.null(z_label)) g = g + xlab(z_label)
  if (!is.null(yrange)) g = g +  scale_y_continuous(limits = yrange)
  print(g)
  return(g)
}
