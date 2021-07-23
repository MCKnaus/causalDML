#' This function implements an ensemble learner with the possibility to extract the
#' smoother matrix of linear smoothers.
#'
#' @param ml List of methods build via \code{\link{create_method}} to be used in ensemble.
#' @param x Covariate matrix of training sample.
#' @param y Vector of outcomes of training sample.
#' @param xnew Covariate matrix of test sample.
#' @param nfolds Number of folds used in cross-validation of ensemble weights (default \code{nfolds=5}).
#' @param weights If TRUE, weights underlying the ensemble prediction for xnew calculated.
#' @param quiet If FALSE, print method that is currently computed.
#'
#' @return Ensemble object containing:
#' \item{ensemble}{the predictions of the trained ensemble for the test sample}
#' \item{best}{the predictions of the best performing single method for the test sample}
#' \item{fit_full}{matrix containing predictions of all considered methods}
#' \item{weights}{If \code{weights=TRUE}, smoother matrix of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver
#' ensemble predictions where each row gives the weight that each training outcome received in the prediction}
#'  \item{nnls_weights}{the weights that each method receives in the ensemble}
#'  \item{mse_cv}{cross-validated MSEs for each method}
#'  \item{fit_cv}{matrix containing the cross-validated predictions for each method}
#'
#'  @export
#'
ensemble = function(ml,
                    x,y,xnew,
                    nfolds=5,
                    weights=FALSE,
                    quiet=TRUE) {
  # Matrix to store the cross-validated predictions
  fit_cv = matrix(NA,nrow(x),length(ml))
  colnames(fit_cv) = sprintf("Method%s",seq(1:length(ml)))
  for (i in 1:length(ml)) {
    if (!is.null(ml[[i]]$name)) colnames(fit_cv)[i] = ml[[i]]$name
  }

  # Get CV folds
  cvf = prep_cf_mat(length(y),nfolds) # see utils_ensemble.R

  # Loop over different folds
  if (length(ml) > 1) {
    for (i in 1:nfolds) {
      # Define training and test sample for this fold
      fold = cvf[,i]
      x_tr = x[!fold,]
      y_tr = y[!fold]
      x_te = x[fold,]

      # Get predictions of all methods for test sample
      fit_cv[fold,] = ensemble_core(ml,x_tr,y_tr,x_te,quiet=quiet)$predictions
    }
    # Cross-validated MSE
    fit_cv[is.na(fit_cv)] = mean(y) # e.g. glmnet produced sometimes NaN for logistic Ridge
    mse_cv = colMeans((c(y) - fit_cv)^2)

    # Calculate weights each method receives
    nnls_weights = nnls::nnls(fit_cv,y)$x
    # In case of perfectly agreeing predictions nnls provides only zeros
    if (sum(nnls_weights) == 0) nnls_weights = nnls_weights + 1 / length(nnls_weights)
    nnls_weights = nnls_weights / sum(nnls_weights)

    # Run all methods on the full sample
    fit_full = ensemble_core(ml,x,y,xnew,weights=weights,quiet=quiet)
    best = fit_full$predictions[,which.min(mse_cv)]
    ensemble = fit_full$predictions %*% nnls_weights

    # Calculate Smoothing matrix if weights=TRUE
    w = NULL
    if (isTRUE(weights)) {
      w = matrix(0,nrow(xnew),nrow(x))
      for (i in 1:length(ml)) {
        w = w + nnls_weights[i] * fit_full$weights[[i]]
      }
      w = Matrix::Matrix(w,sparse=T)
    }

    colnames(fit_full$predictions) = names(mse_cv) = names(nnls_weights) = colnames(fit_cv)
  }
  else { # in case only one method defined, no cross-validation needed
    fit_full = ensemble_core(ml,x,y,xnew,weights=weights,quiet=quiet)
    ensemble = best = fit_full$predictions
    w = nnls_weights = mse_cv = fit_cv = NULL
    if (isTRUE(weights)) w = fit_full$weights[[1]]
  }

  output = list("ensemble" = ensemble,"best" = best,"fit_full" = fit_full,"weights" = w,
       "nnls_weights" = nnls_weights, "mse_cv" = mse_cv, "fit_cv" = fit_cv)
  class(output) = "ensemble"
  output
}


#' Core function of \code{\link{ensemble}}.
#'
#' @param ml List of methods build via \code{\link{create_method}} to be used in ensemble.
#' @param x_tr Covariate matrix of training sample.
#' @param y_tr Vector of outcomes of training sample.
#' @param x_te Covariate matrix of test sample.
#' @param weights If TRUE, weights underlying the ensemble prediction for x_te calculated.
#' @param quiet If FALSE, print method that is currently computed.
#'
#' @keywords internal
#'
#' @return Returns list containing:
#' \item{predictions}{\code{nrow(x_te)} x \code{length(ml)} matrix with the predictions of each method}
#' \item{weights}{If \code{weights=TRUE}, list of length(ml) with smoother matrices of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions of each method where each row gives the weight that each training
#' outcome received in the prediction}
#'
ensemble_core = function(ml,
                         x_tr,y_tr,x_te,
                         weights=FALSE,
                         quiet=TRUE) {
  # Initialize objects to be filled
  fit_mat = matrix(NA,nrow(x_te),length(ml))
  weights_list = vector("list",length(ml))

  # Loop over specified methods
  for (i in 1:length(ml)) {
    wrapper = paste0(ml[[i]]$method,"_fit")
    if (isFALSE(quiet)) print(wrapper)

    # Check whether subset of variables specified and/or additional arguments are defined and run method
    if (is.null(ml[[i]]$x_select) & length(ml[[i]]$args) == 0)          fit = do.call(wrapper,list(x=x_tr,y=y_tr))
    else if (is.null(ml[[i]]$x_select) & !(length(ml[[i]]$args) == 0))  fit = do.call(wrapper,list(x=x_tr,y=y_tr,args=ml[[i]]$args))
    else if (!is.null(ml[[i]]$x_select) & length(ml[[i]]$args) == 0)    fit = do.call(wrapper,list(x=x_tr[,ml[[i]]$x_select],y=y_tr))
    else                                                                fit = do.call(wrapper,list(x=x_tr[,ml[[i]]$x_select],y=y_tr,args=ml[[i]]$args))
    # Get predictions
    if (is.null(ml[[i]]$x_select))  temp = do.call(paste0("predict.",wrapper),list(fit,x=x_tr,y=y_tr,xnew=x_te,weights=weights))
    else                            temp = do.call(paste0("predict.",wrapper),list(fit,x=x_tr[,ml[[i]]$x_select],y=y_tr,xnew=x_te[,ml[[i]]$x_select],weights=weights))

    fit_mat[,i] = temp$prediction
    weights_list[[i]] = temp$weights
  }

  list("predictions" = fit_mat, "weights" = weights_list)
}


#' Creates the methods to be used in \code{\link{ensemble}}
#'
#' @param method Choose method from currently \code{c("mean","ols",ridge","plasso",forest_grf","lasso")}.
#' @param x_select Optional logical vector of length equal to the number of columns of the covariate matrix
#' indicating which variables should be used by this method. E.g. tree-based methods usually should not be provided
#' with the interactions that Lasso is using.
#' @param args Optional list containing the additional arguments that should be passed to the underlying method.
#' @param name Optional string naming the method.
#'
#' @return Method that can be passed to \code{\link{ensemble}}.
#'
#' @export
#'
create_method = function(method,
                         x_select=NULL,
                         args=list(),
                         name=NULL) {

  if (!(is.character(method) & length(method) == 1)) stop("Provide single string to define method.")
  if (!(any(method == c("mean","ols","ridge","plasso","forest_grf","lasso")))
  ) stop("Provide one of these options c(\"mean\",\"ols\",\"ridge\",\"plasso\",\"forest_grf\",\"lasso\").")
  if (!(is.null(args) | is.list(args))) stop("Provide either NULL or list for args.")
  if (!(is.null(x_select) | is.logical(x_select))) stop("Provide either NULL or logical for x_select.")
  if (!((is.character(name) & length(name) == 1) | is.null(name))) stop("Provide single string to name method.")

  list(method=method,args=args,x_select=x_select,name=name,weights=weights)
}
