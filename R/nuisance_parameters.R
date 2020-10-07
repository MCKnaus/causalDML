#' Cross-fitted ensemble prediction of propensity score nuisance parameter.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param x Covariate matrix.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param cv Number of cross-validation when estimating ensemble (default 5).
#' @param path Optional path to save the \code{\link{ensemble}} objects for later inspection.
#' Saved as Ensemble_Wi where i is the number of the treatment in multiple treatment settings.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return Returns n x T+1 matrix with each columns containing the propensity score for the treatment corresponding to w_mat.
#'
#' @export
#'
nuisance_e = function(ml,w_mat,x,cf_mat,
                    cv=5,
                    path = NULL,
                    quiet=TRUE) {
  if (isFALSE(quiet)) print("Propensity score")

  # Initialize nuisance matrix
  e_mat = matrix(NA,nrow(w_mat),ncol(w_mat))
  colnames(e_mat) = colnames(w_mat)

  ## Calculate treatment probabilities
  # Binary treatment
  if (ncol(w_mat) == 2) {
    if (!is.null(path)) path_tem = paste0(path,"Ensemble_W")
    else path_tem = NULL
    e_mat[,1] = nuisance_cf(ml,w_mat[,1],x,cf_mat,cv=cv,path=path_tem,quiet=quiet)
    e_mat[,2] = 1 - e_mat[,1]
  }
  # Multiple treatment
  else if (ncol(w_mat) > 2) {
    for (i in 1:ncol(w_mat)) {
      if (!is.null(path)) path_tem = paste0(path,"Ensemble_W",toString(i))
      else path_tem = NULL
      e_mat[,i] = nuisance_cf(ml,w_mat[,i],x,cf_mat,cv=cv,path=path_tem,quiet=quiet)
    }
    e_mat = e_mat / rowSums(e_mat)
  }
  else stop("Provide treatment indicator matrix with at least 2 columns")

  return(e_mat)
}


#' Cross-fitted ensemble prediction of outcome regression nuisance parameter.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}.
#' @param y Numerical vector containing the outcome variable.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param x Covariate matrix.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param cv Number of cross-validation when estimating ensemble (default 5).
#' @param weights If TRUE, prediction weights of the outcome nuisance extracted and saved (requires to provide a path).
#' @param path Optional path to save the \code{\link{ensemble}} objects for later processing.
#' Saved as Ensemble_Yi where i is the number of the treatment in multiple treatment settings.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return Returns n x T+1 matrix with each columns containing the predicted outcome for the treatment corresponding to w_mat.
#'
#' @export
#'
nuisance_m = function(ml,y,w_mat,x,cf_mat,
                      cv=5,
                      weights=FALSE,
                      path=NULL,
                      quiet=TRUE) {

  if (isFALSE(quiet)) print("Outcome regression")

  # Initialize nuisance matrix
  m_mat = matrix(NA,nrow(w_mat),ncol(w_mat))
  colnames(m_mat) = colnames(w_mat)

  ## Calculate outcome predictions
  for (i in 1:ncol(w_mat)) {
    if (!is.null(path)) path_tem = paste0(path,"Ensemble_Y",toString(i))
    else path_tem = NULL
    m_mat[,i] = nuisance_cf(ml,y,x,cf_mat,cv=cv,subset=(w_mat[,i]),weights=weights,path=path_tem,quiet=quiet)
  }

  return(m_mat)
}



#' Cross-fitting of nuisance parameter with \code{\link{ensemble}}.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation.
#' Methods can be created by \code{\link{create_method}}.
#' @param y Vector of variable to be predicted.
#' @param x Matrix of covariates.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param cv Number of cross-validation when estimating ensemble (default 5).
#' @param subset Optional logical vector if only subset of data should be used for prediction.
#' @param weights If TRUE, prediction weights of the outcome nuisance extracted and saved (requires to provide a path).
#' @param path Optional path to save the \code{\link{ensemble}} of each fold for later processing.
#' Saved as path + "_foldi" where i is the fold number.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return Returns a n x 1 matrix of nuisance parameters.
#'
#' @keywords internal
#'
nuisance_cf = function(ml,y,x,cf_mat,
                       cv=5,
                       subset=NULL,
                       weights=FALSE,
                       path=NULL,
                       quiet=TRUE) {
  # Checks
  if (is.numeric(cf_mat)) {
    if (!all(cf_mat %in% 0:1) | !is.matrix(cf_mat)) stop("Please provide cf_mat as binary indicator matrix. E.g. use function prep_cf_mat")
    if (nrow(cf_mat) != length(y)) stop("cf_mat indicator matrix nrows different from # of obs.")
    if (length(y) != sum(cf_mat)) stop("cf_mat indicator matrix does not sum to number of observations.")
  }
  if (isTRUE(weights) & is.null(path)) stop("Provide path if weights=TRUE to save ensemble objects with weights for
                                            later processing or set weights=FALSE.")

  if (is.null(subset)) subset = rep(TRUE,length(y))

  np = rep(NA,length(y))

  for (i in 1:ncol(cf_mat)) {
    if (isFALSE(quiet)) print(paste("Cross-fitting fold:",toString(i)))
    fold = cf_mat[,i]
    x_tr = x[!fold & subset,]
    y_tr = y[!fold & subset]
    x_te = x[fold,]

    ens = do.call(ensemble,c(list(ml=ml,x=x_tr,y=y_tr,xnew=x_te,nfolds=cv,weights=weights,quiet=quiet)))
    np[fold] = ens$ensemble
    if (!is.null(path)) {
      save(ens, file = paste0(path,"_fold",toString(i),".RData"))
    }
    rm(ens); gc()
  }
  return(np)
}
