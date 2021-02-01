#' This function uses the \code{\link{glmnet}} package to estimate the coefficient paths
#' and cross-validates least squares Lasso AND Post-Lasso
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param w vector of weights
#' @param kf number of folds in k-fold CV
#' @param ... Pass \code{\link{glmnet}} options
#' @import glmnet
#'
#' @return List with the names of selected variables at cross-validated minima for Lasso and Post-Lasso
#'
plasso = function(x,y,
                  w=NULL,
                  kf = 10,
                  ...) {

  # Handle potentially provided sample weights, otherwise create weight vector of ones
  w = handle_weights(w,nrow(x))

  # Lasso with full estimation sample
  lasso_full = glmnet(x,y,weights = as.vector(w),family="gaussian",...)
  coef_lasso_full = coef(lasso_full)                   # Save coefficients to extract later the ones at the CV minima
  nm_coef_lasso_full = rownames(coef_lasso_full)[-1]   # Save variable names that were used and kick out intercept
  lambda = lasso_full$lambda                           # Save grid to use the same in cross-validation

  ###################################################
  ### Cross-validation with Lasso and Post Lasso ####
  ###################################################

  ## Figure out variables that were always inactive to drop them in CV
  # Figure out which coefficients are inactive at each Lasso grid
  inact_coef_full = (coef_lasso_full == 0)         # boolean matrix

  # Initialize matrices for MSE of Lasso and post-lasso for each grid point and CV fold
  cv_MSE_lasso = matrix(nrow = kf,ncol = length(lambda))
  cv_MSE_plasso = matrix(nrow = kf,ncol = length(lambda))

  # Get indicator for CV samples
  split = stats::runif(nrow(x))
  cvgroup = as.numeric(cut(split,stats::quantile(split, probs = seq(0,1,1/kf)),include.lowest = TRUE))  # groups for K-fold cv
  list = 1:kf                            # Needed later in the loop to get the appropriate samples


  ## Start loop for cross-validation of Lasso and Post-Lasso
  for (i in 1:kf) {
    CV = CV_core(x,y,w,cvgroup,list,i,lambda,...)

    ## Extract MSE of Lasso and Post-Lasso
    cv_MSE_lasso[i,] = CV$MSE_lasso
    cv_MSE_plasso[i,] = CV$MSE_plasso
  }       # end loop over folds

  ## Calculate mean MSEs over all folds
  cv_MSE_lasso[is.na(cv_MSE_lasso)] = max(cv_MSE_lasso) # can happen if glmnet does not go over the full grid
  cv_MSE_plasso[is.na(cv_MSE_plasso)] = max(cv_MSE_plasso) # and/or Post-Lasso has not full rank
  mean_MSE_lasso = colMeans(cv_MSE_lasso)
  mean_MSE_plasso = colMeans(cv_MSE_plasso)
  mean_MSE_lasso[is.na(mean_MSE_lasso)] = max(mean_MSE_lasso,na.rm=T) + 1e-7 # can happen if glmnet does not go over the full grid
  mean_MSE_plasso[is.na(mean_MSE_plasso)] = max(mean_MSE_plasso,na.rm=T) + 1e-7 # and/or Post-Lasso has not full rank

  ## Get grid position of minimum MSE
  ind_min_l = which.min(mean_MSE_lasso)
  ind_min_pl = which.min(mean_MSE_plasso)
  # Get names at minima
  names_l = names(coef_lasso_full[,ind_min_l])[which(coef_lasso_full[,ind_min_l] != 0)]
  names_pl = names(coef_lasso_full[,ind_min_pl])[which(coef_lasso_full[,ind_min_pl] != 0)]

  ## Get Lasso coefficients at minimum of Lasso
  coef_min_l = coef_lasso_full[,ind_min_l][which(coef_lasso_full[,ind_min_l] != 0)]
  ## Get Lasso coefficients at minimum of Post-Lasso
  coef_min_pl = coef_lasso_full[,ind_min_pl][which(coef_lasso_full[,ind_min_pl] != 0)]

  ## Return names and coefficients
  output = list("lasso_full"=lasso_full,"kf"=kf,
                "cv_MSE_lasso"=cv_MSE_lasso,"cv_MSE_plasso"=cv_MSE_plasso,
                "mean_MSE_lasso" = mean_MSE_lasso, "mean_MSE_plasso" = mean_MSE_plasso,
                "ind_min_l" = ind_min_l,"ind_min_pl" = ind_min_pl,
                "lambda_min_l" = lambda[ind_min_l],"lambda_min_pl" = lambda[ind_min_pl],
                "names_l" = names_l,"names_pl" = names_pl,
                "coef_min_l" = coef_min_l,"coef_min_pl" = coef_min_pl)

  class(output) = "plasso"
  output
}


#' Predict after Post-Lasso.
#'
#' @param plasso \code{\link{plasso}} object
#' @param x Covariate matrix that was used for training
#' @param y Outcomes that were used for training
#' @param xnew Matrix of new values for x at which predictions are to be made
#' @param se_rule If equal to zero predictions from CV minimum (default). Negative values go in the direction of smaller
#' models (e.g. se_rule=-1 creates the standard 1SE rule), positive values go in the direction of larger models
#' (e.g. se_rule=1 creates the standard 1SE+ rule)
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
predict.plasso = function(plasso,x,y,
                          xnew=NULL,
                          se_rule=0,
                          weights=FALSE) {

  if (is.null(xnew)) xnew = x
  x = add_intercept(x)
  xnew = add_intercept(xnew)

  # Standard error of folds
  oneSE_lasso = sqrt(apply(plasso$cv_MSE_lasso, 2, var)/plasso$kf)
  oneSE_plasso = sqrt(apply(plasso$cv_MSE_plasso, 2, var)/plasso$kf)
  oneSE_lasso[is.na(oneSE_lasso)] = 0
  oneSE_plasso[is.na(oneSE_plasso)] = 0

  # Find Lambda
  ind_Xse_l = find_Xse_ind(plasso$mean_MSE_lasso,plasso$ind_min_l,oneSE_lasso,se_rule)
  ind_Xse_pl = find_Xse_ind(plasso$mean_MSE_plasso,plasso$ind_min_pl,oneSE_plasso,se_rule)

  # Fitted values for lasso
  fit_lasso = xnew %*% coef(plasso$lasso_full)[,ind_Xse_l]

  # Fitted values for post lasso
  nm_act = names(coef(plasso$lasso_full)[,ind_Xse_pl])[which(coef(plasso$lasso_full)[,ind_Xse_pl] != 0)]
  xact = x[,nm_act]
  xactnew = xnew[,nm_act]
  hat_mat = xactnew %*% solve(crossprod(xact)) %*% t(xact)
  fit_plasso = hat_mat %*% y
  if (weights==FALSE) hat_mat = NULL

  list("lasso"=fit_lasso,"plasso"=fit_plasso,"weights"=hat_mat)
}


#' Summary of Post-Lasso model
#'
#' @param plasso \code{\link{plasso}} object
#'
#' @return Prints cross-validated MSE and active variables for Lasso and Post-Lasso.
#'
summary.plasso = function(plasso) {
  # Comparison of minimum MSE
  cat("\n\n Minimum CV MSE Lasso:",toString(min(plasso$mean_MSE_lasso,na.rm = TRUE)))
  cat("\n\n Minimum CV MSE Post-Lasso:",toString(min(plasso$mean_MSE_plasso,na.rm = TRUE)))

  # Show names of active variables at respective minima
  cat("\n\n Active variables at CV minimum of Lasso: \n")
  print(names(coef(plasso$lasso_full)[,plasso$ind_min_l])[which(coef(plasso$lasso_full)[,plasso$ind_min_l] != 0)])

  cat("\n\n Active variables at CV minimum of Post-Lasso: \n")
  print(names(coef(plasso$lasso_full)[,plasso$ind_min_pl])[which(coef(plasso$lasso_full)[,plasso$ind_min_pl] != 0)])
}


#' Plot of cross-validation curves.
#'
#' @param plasso \code{\link{plasso}} object
#'
plot.plasso = function(plasso) {

  # Standard error of folds
  oneSE_lasso = sqrt(apply(plasso$cv_MSE_lasso, 2, var)/plasso$kf)
  oneSE_plasso = sqrt(apply(plasso$cv_MSE_plasso, 2, var)/plasso$kf)
  oneSE_lasso[is.na(oneSE_lasso)] = 0
  oneSE_plasso[is.na(oneSE_plasso)] = 0

  # Calculate 1SE bands
  lasso_1se_up = plasso$mean_MSE_lasso+oneSE_lasso
  lasso_1se_low = plasso$mean_MSE_lasso-oneSE_lasso
  plasso_1se_up = plasso$mean_MSE_plasso+oneSE_plasso
  plasso_1se_low = plasso$mean_MSE_plasso-oneSE_plasso

  # Get ranges for the graph
  xrange = range(log(plasso$lasso_full$lambda))
  yrange = c(max(-1.7e+308,min(lasso_1se_low,plasso_1se_low)),
             min(1.7e+308,max(lasso_1se_up,plasso_1se_up)))

  # Plot mean lines
  ylab = "Mean-squared Error"
  graphics::plot(xrange,yrange,type="n",xlab="Log Lambda",ylab=ylab)
  graphics::lines(log(plasso$lasso_full$lambda),plasso$mean_MSE_lasso,lwd=1.5,col="blue")
  graphics::lines(log(plasso$lasso_full$lambda),plasso$mean_MSE_plasso,lwd=1.5,col="red")

  # Plot upper and lower 1SE lines
  graphics::lines(log(plasso$lasso_full$lambda),lasso_1se_up,lty=2,lwd=1,col="blue")
  graphics::lines(log(plasso$lasso_full$lambda),lasso_1se_low,lty=2,lwd=1,col="blue")
  graphics::lines(log(plasso$lasso_full$lambda),plasso_1se_up,lty=2,lwd=1,col="red")
  graphics::lines(log(plasso$lasso_full$lambda),plasso_1se_low,lty=2,lwd=1,col="red")

  # Show location of minima
  graphics::abline(v = log(plasso$lambda_min_l), lty = 1, col="blue")
  graphics::abline(v = log(plasso$lambda_min_pl), lty = 1, col="red")

  # Print legend
  graphics::legend('top',c("Lasso MSE","Lasso MSE+-1SE","Post-Lasso MSE","Post-Lasso MSE+-1SE","# active coeff."), lty = c(1,2,1,2,1),
                   col=c('blue','blue','red','red','forestgreen'),ncol=1,bty ="n",cex=0.7)

  # Open a new graph for number of coefficients to be written on existing
  graphics::par(new = TRUE)
  graphics::plot(log(plasso$lasso_full$lambda),plasso$lasso_full$df, axes=F, xlab=NA, ylab=NA, cex=1.2,type="l",col="forestgreen",lwd=1.5)
  graphics::axis(side = 4)
  graphics::mtext(side = 4, line = 3, "# active coefficients")
}


#' This function contains the core parts of the CV for Lasso and Post-Lasso

#' @param x covariate matrix to be used in CV
#' @param y vector of outcomes
#' @param w vector of weight
#' @param cvgroup categorical with k groups to identify folds
#' @param list list 1:k
#' @param i number of fold that is used for prediction
#' @param lambda series of lambdas used
#' @param ... Pass \code{\link{glmnet}} options
#'
#' @return MSE_lasso / MSE_plasso: means squared errors for each lambda
#'
#' @keywords internal
#'

CV_core = function(x,y,w,cvgroup,list,i,lambda,...) {

  # Get estimation and prediction sample for this specific fold
  x_est_cv = subset(x,cvgroup %in% list[-i])
  y_est_cv = subset(y,cvgroup %in% list[-i])
  w_est_cv = subset(w,cvgroup %in% list[-i])
  x_pred_cv = subset(x,cvgroup %in% list[i])
  y_pred_cv = subset(y,cvgroup %in% list[i])
  w_pred_cv = subset(w,cvgroup %in% list[i])

  # Normalize the weights to N
  w_est_cv = norm_w_to_n(w_est_cv)
  w_pred_cv = norm_w_to_n(w_pred_cv)

  # Estimate Lasso for this fold using the grid of the full sample
  lasso_cv = glmnet(x_est_cv, y_est_cv,lambda = lambda,weights=as.vector(w_est_cv),
                    family="gaussian",...)
  coef_lasso_cv = coef(lasso_cv)                                       # Save coefficients at each grid point

  # Predicted values with lasso coefficients for prediction sample at each grid
  fit_lasso = predict(lasso_cv,x_pred_cv)
  if (ncol(fit_lasso) != length(lambda)) {
    fit_lasso = cbind(fit_lasso,matrix(NA,nrow(fit_lasso),(length(lambda)-ncol(fit_lasso))))
  }

  ## Now calculate predicted values for post Lasso at each grid
  # Initialize first matrix for fitted values
  fit_plasso = matrix(NA,nrow = nrow(fit_lasso),ncol = ncol(fit_lasso))

  # Figure out which coefficients are active at each Lasso grid
  act_coef = (coef_lasso_cv != 0)         # boolean matrix

  ## Calculate full covariance matrix X'X and X'y once and select only the relevant in the loop below (much faster)
  # First, figure out variables that were active at least once and get only these
  act_once = apply(act_coef,1,any)
  nm_all_act_coef = rownames(act_coef)[act_once]
  if (length(nm_all_act_coef) == 1) {
    warning("No variables selected in one CV, even for lowest lambda, might reconsider choice of penalty term grid")
    fit_plasso = fit_lasso
  } else {
    ## Create the covariate matrix to "manually" calculate the fitted values (much faster than the build in lm command)
    if (sum(abs(coef_lasso_cv[1,]))==0) {   # Indicates that no intercept was used
      x_all_act = x_est_cv[,nm_all_act_coef]
    }
    else if (sum(abs(coef_lasso_cv[1,]))!=0) {    # Indicates that intercept was used
      x_all_act = add_intercept(x_est_cv[,nm_all_act_coef[2:length(nm_all_act_coef)],drop=F])
      colnames(x_all_act)[1] = "(Intercept)"
      # add intercept also to prediction sample
      x_pred_cv = add_intercept(x_pred_cv[,nm_all_act_coef[2:length(nm_all_act_coef)],drop=F])
      colnames(x_pred_cv)[1] = "(Intercept)"
    }
    else {
      stop("Something strange happens with the intercepts at the Lasso path")
    }
    # Add weights
    x_w = apply(x_all_act,2,`*`,sqrt(w_est_cv))
    y_w = y_est_cv * sqrt(w_est_cv)

    # Get X'X
    XtX_all = crossprod(x_w)
    # Get X'y
    Xty_all = crossprod(x_w,y_w)

    # Initialize vector to save chosen variable names to figure out whether new variables were added from the last to the next grid
    nm_act_coef_prev = "nobody should call a variable like this"

    ## Loop over the grid of Lambdas to get the Post-Lasso predictions
    for (j in 1:length(lambda)) {
      # Get names of active variables at this grid
      nm_act_coef = rownames(act_coef)[act_coef[,j]]
      if (identical(nm_act_coef,character(0))==TRUE) {
        # print("No variable selected at this grid => no prediction")
        next
      }
      if (identical(nm_act_coef,nm_act_coef_prev) == TRUE & j>1) {
        # print("Same variables selected as in previous grid => Post-Lasso predictions remain unchanged")
        fit_plasso[,j] = fit_plasso[,(j-1)]
        next
      }
      else {     # Get prediction covariate matrix for that grid
        x_ols_pred = as.matrix(x_pred_cv[,nm_act_coef])
      }

      # Get OLS Post-Lasso predictions for this grid point
      fit = fitted_values(XtX_all,Xty_all,x_ols_pred,nm_act_coef)
      if (is.null(fit) & j==1) {
        fit_plasso[,j] = rep(mean(y == 1),nrow(fit_plasso))
        next
      }
      if (is.null(fit) & j>1) {
        # cat("\n X'X not invertible at grid",toString(j),": Use last feasible coefficients")
        fit_plasso[,j] = fit_plasso[,(j-1)]
        next
      } else{
        fit_plasso[,j] = fit
      }

      # Update current active covariates for the check whether it changes in the next grid
      nm_act_coef_prev = nm_act_coef
    }       # end loop over grids

  } # end if at least one var selected

  # Matrix with "real" outcome values for each grid
  y_rep = matrix(rep(y_pred_cv,length(lambda)),nrow = nrow(fit_lasso),ncol = ncol(fit_lasso))

  # Get RMSE
  SE_lasso = (y_rep - fit_lasso)^2
  SE_plasso = (y_rep - fit_plasso)^2
  if (!is.null(w)) {                                     # Weight errors by "sampling weights"
    SE_lasso = apply(SE_lasso,2,"*",w_pred_cv)
    SE_plasso = apply(SE_plasso,2,"*",w_pred_cv)
  }
  MSE_lasso = apply(SE_lasso,2,mean)
  MSE_plasso = apply(SE_plasso,2,mean)

  list("MSE_lasso" = MSE_lasso,"MSE_plasso" = MSE_plasso)
}




#' This helper function extracts a subset of active variables (nm_act) of the
#' relavant variables from X'X and X'y to get out-of-sample predictions
#' for a matrix containing only the active variables
#' This speeds up the cross-validation for post-Lasso to a large extent
#'
#' @param XtX_all Crossproduct of all covariates
#' @param Xty_all Crossproduct of covariates and outcome
#' @param x_pred Covariates matrix of the PREDICTION sample
#' @param nm_act names of active variables
#'
#' @return Fitted values in the prediction sample
#'
#' @keywords internal
#'

fitted_values = function (XtX_all,Xty_all,x_pred,nm_act) {

  # Extract relevant rows and columns
  XtX = XtX_all[nm_act,nm_act]
  Xty = Xty_all[nm_act,]

  # Calculate point estimates in estimation sample
  b = tryCatch(solve(XtX, Xty), error=function(e) NULL) # much faster than standard solve
  if (is.null(b)) {        # In case b has not been estimated
    fit_val = b
  } else {
    # Calculate fitted values in prediction sample
    fit_val = x_pred%*%b
  }
  return(fit_val)
}


#' Function to normalize weights to N or to N in treated and controls separately
#' @param w vector of weights that should be normalized
#' @param d vector of treament indicators
#'
#' @return Normalized weights
#'
#' @keywords internal
#'

norm_w_to_n = function(w,d=NULL) {

  if (is.null(d)) {
    ## Normalize weights to sum up to N
    w = w / sum(w)* nrow(w)
  } else {
    # Separate weights of treated and controls
    w1 = w * d
    w0 = w * (1-d)
    # Normalize weights to sum to N in both groups
    w1 = w1 / sum(w1) * nrow(w)
    w0 = w0 / sum(w0) * nrow(w)
    # Unify weights again
    w = w1 + w0
  }
  return(w)
}


#' Helper function finds the position for pre-specified SE rules
#'
#' @param CV Vector of cross-validated criterion
#' @param ind_min Index of cross-validated minimum
#' @param oneSE Standard error of cross-validated criterion at the minimum
#' @param factor Factor in which direction to go. Negative smaller model, positive larger model
#'
#' @return Index on the Lambda grid
#'
#' @keywords internal
#'
find_Xse_ind = function(CV,ind_min,oneSE,factor) {
  cv_temp = CV - (CV[ind_min] + abs(factor) * oneSE[ind_min])
  if (factor < 0) {
    for (i in ind_min:1) {
      ind = i
      if (cv_temp[i] < 0) next
      else if (cv_temp[i] > 0) break
    }
  } else if (factor < 0) {
    for (i in ind_min:length(oneSE)) {
      ind = i
      if (cv_temp[i] < 0) next
      else if (cv_temp[i] > 0) break
    }
  } else ind = ind_min
  return(ind)
}


#' Adds an intercept to a matrix
#' @param mat Any matrix
#'
#' @return Matrix with intercept
#' @keywords internal
#'
add_intercept = function(mat) {
  if (is.null(dim(mat))) mat = as.matrix(mat,ncol=1)
  mat = cbind(rep(1,nrow(mat)),mat)
  colnames(mat)[1] = "(Intercept)"
  return(mat)
}


#' Sanitizes potential sample weights
#' @param w Vector of weights or null if no weights provided
#' @param n number of observations
#'
#' @return Vector of weights
#'
#' @keywords internal
#'
handle_weights = function(w,n) {
  # Create weights of ones if no weights are specified
  if (is.null(w)) {
    w <- as.matrix(rep(1,n),nrow = n,ncol = 1)
  } else {
    w <- as.matrix(w,nrow = n,ncol = 1)
  }
  colnames(w) <- "w"

  # Normalize the weights either to N
  w <- norm_w_to_n(w)
  return(w)
}

