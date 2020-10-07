#' Creates matrix of binary treament indicators (n x T+1) where each column contains
#' a binary indicator for one treatment.
#'
#' @param w Treatment vector. Provide as factor to control ordering of the treatments,
#' otherwise program orders treatments in ascending order or alphabetically.
#'
#' @return Logical matrix of treatment indicators (n x T+1).
#'
#' @import stats
#'
#' @export
#'
prep_w_mat = function(w) {
  # Checks
  if (length(unique(w)) <= 1) stop("Need at least two values in treatment vector")
  if (!is.factor(w)) w = factor(w, levels = sort(unique(w)))

  # Create one-hot matrix for each category
  w_mat = model.matrix(~0+w)
  colnames(w_mat) = gsub("w","",colnames(w_mat))
  return(w_mat == 1)
}


#' Creates matrix of binary cross-fitting fold indicators (n x # cross-folds)
#'
#' @param n Number of observations.
#' @param cf Number of cross-fitting folds.
#' @param cl Optional vector of cluster variable if cross-fitting should account for clusters.
#'
#' @importFrom dplyr ntile
#' @import stats
#'
#' @return Logical matrix of cross-fitting folds (n x # folds).
#'
#' @export
#'
prep_cf_mat = function(n,cf,cl=NULL) {

  if (cf == 1) cf_mat = matrix(rep(1,n),ncol=1)

  if (!is.null(cl)) {
    fold = ntile(runif(length(unique(cl))),cf)
    fold = factor(fold[match(cl,unique(cl))])
    cf_mat = model.matrix(~0+fold)
  }
  else {
    fold = factor(ntile(runif(n),cf))
    cf_mat = model.matrix(~0+fold)
  }
  colnames(cf_mat) = sprintf("CF %d",1:cf)

  return(cf_mat == 1)
}


#' Runs classification analysis (CLAN) for CATEs.
#'
#' @param cate Vector of CATEs.
#' @param x Covariate matrix.
#' @param q Number of splits (defaut 5).
#'
#' @importFrom dplyr ntile
#' @import matrixStats
#'
#' @return Matrix containing means of covariate in each group.
#'
#' @export
#'
clan = function(cate,x,q=5) {
  qt = ntile(cate,q)
  x = cbind(cate,x)
  res = matrix(NA,ncol(x),q)
  rownames(res) = colnames(x)
  for (i in 1:q) {
    res[,i] = colMeans(x[qt==i,])
  }
  return(res)
}
