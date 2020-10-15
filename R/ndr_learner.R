#' DR- and NDR-learner.
#'
#' This function produces out-of-sample predictions of conditional average treatment effects (CATEs) for
#' all individuals in the sample using the DR-learner (Kennedy, 2020) and the normalized DR-learner (Knaus, 2020).
#' The involved predictions are based on \code{\link{ensemble}}. See algorithms 1 and 2 in Knaus (2020) for details.
#'
#' @param y Numerical vector containing the outcome variable.
#' @param w Treatment vector. Provide as factor to control ordering of the treatments,
#' otherwise program orders treatments in ascending order or alphabetically.
#' @param x Covariate matrix.
#' @param ml_w List of methods to be used in ensemble estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param ml_y List of methods to be used in ensemble estimation of outcome regression.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param ml_tau List of methods to be used in ensemble estimation of CATEs.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param compare_all Relevant multiple treatments: If FALSE, only comparisons to first treatment considered.
#' @param nfolds Number of folds used in cross-validation of ensemble weights (default \code{nfolds=5}).
#' @param path Optional path to save the \code{\link{ensemble}} objects used to predict nuisance parameters and CATEs
#' for later processing. Outputs of each fold i saved in new subfolder NDR_foldi.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return  \item{cates}{Returns n x 2 matrix containing DR- and NDR-learner predictions in case of a binary treatment or
#' a number of comparisons x n x 2 array with the predictions.}
#'          \item{list}{A list of the four \code{\link{ndr_oos}} outputs.}
#'
#' @references
#' \itemize{
#' \item Kennedy, E. H. (2020). Optimal doubly robust estimation of heterogeneous causal effects.
#' arXiv preprint arXiv:2004.14497. \url{http://arxiv.org/abs/2004.14497}
#' \item Knaus, M. C. (2020). Double machine learning based program evaluation under unconfoundedness.
#'   arXiv preprint arXiv:2003.03191.\url{http://arxiv.org/abs/2003.03191}
#' }
#'
#' @export
#'
ndr_learner = function(y,w,x,
                       ml_w = list(create_method("forest_grf")),
                       ml_y = list(create_method("forest_grf")),
                       ml_tau = list(create_method("forest_grf")),
                       compare_all = TRUE,
                       nfolds=5,
                       path=NULL,
                       quiet=TRUE) {
  # Create indicator matrix
  wm = prep_w_mat(w)

  # Split sample in four folds
  cfm = prep_cf_mat(length(y),4)

  # How many comparisons are specified?
  if (isTRUE(compare_all)) num_comp = ncol(wm)*(ncol(wm)-1)/2
  else num_comp = ncol(wm)-1

  # Initialize vectors to store results
  list_ndr_oos = vector("list",4)
  cates = array(NA,c(num_comp,nrow(x),2))

  # Loop over four folds
  for (i in 1:4) {
    if (isFALSE(quiet)) print(paste("NDR-learner fold", toString(i)))
    oos = cfm[,i]
    if (!is.null(path)) {
      dir.create(paste0(path,"/NDR_fold",toString(i)))
      path_temp = paste0(path,"/NDR_fold",toString(i),"/")
    }
    list_ndr_oos[[i]] = ndr_oos(y[!oos],w[!oos],x[!oos,],ml_w = ml_w,
                                        ml_y = ml_y,ml_tau = ml_tau,
                                        cf_mat = cfm[!oos,-i],compare_all = FALSE,xnew=x[oos,],
                                        path=path_temp,quiet=FALSE)
    for (j in 1:num_comp) {
      cates[j,oos,] = list_ndr_oos[[i]]$cates[[j]]
    }
  }
  if (ncol(wm) == 1) cates = cates[[1]]

  list("cates"=cates,"list"=list_ndr_oos)
}


#' This function produces out-of-sample predictions of conditional average treatment effects (CATEs) using the
#' DR-learner (Kennedy, 2020) and the normalized DR-learner (Knaus, 2020).
#' It executes steps 1 to 4 of algorithms 1 and 2 in Knaus (2020) as part of the \code{\link{ndr_learner}}.
#'
#' @param y Numerical vector containing the outcome variable.
#' @param w Treatment vector. Provide as factor to control ordering of the treatments,
#' otherwise program orders treatments in ascending order or alphabetically.
#' @param x Covariate matrix.
#' @param xnew Covariate matrix of test sample.
#' @param ml_w List of methods to be used in ensemble estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param ml_y List of methods to be used in ensemble estimation of outcome regression.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param ml_tau List of methods to be used in ensemble estimation of CATEs.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param cf_mat Optional logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}). Otherwise created internally.
#' @param compare_all Relevant multiple treatments: If FALSE, only comparisons to first treatment considered.
#' @param nfolds Number of folds used in cross-validation of ensemble weights (default \code{nfolds=5}).
#' @param path Optional path to save the \code{\link{ensemble}} objects for later processing.
#' IATE objects saved in new subfolder Comparisonij where i is the index of the control and j of the treated group.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return  \item{cates}{n x 2 matrix containing DR- and NDR-learner predictions in case of a binary treatment or
#' a list of n x 2 matrices with the specified comparisons.}
#'          \item{APO}{\code{\link{APO_dml}} object containing the underlying nuisance parameters etc.}
#'          \item{ATE}{\code{\link{ATE_dml}} object containing the underlying doubly robust score.}
#'
#' @references
#' \itemize{
#' \item Kennedy, E. H. (2020). Optimal doubly robust estimation of heterogeneous causal effects.
#' arXiv preprint arXiv:2004.14497. \url{http://arxiv.org/abs/2004.14497}
#' \item Knaus, M. C. (2020). Double machine learning based program evaluation under unconfoundedness.
#'   arXiv preprint arXiv:2003.03191.\url{http://arxiv.org/abs/2003.03191}
#' }
#'
#' @export
#'
ndr_oos = function(y,w,x,xnew,
                       ml_w = list(create_method("forest_grf")),
                       ml_y = list(create_method("forest_grf")),
                       ml_tau = list(create_method("forest_grf")),
                       cf_mat = NULL,
                       compare_all = TRUE,
                       nfolds=5,
                       path=NULL,
                       quiet=TRUE) {
  # Create indicator matrices
  if (is.null(cf_mat)) cf_mat = prep_cf_mat(length(y),3)
  wm = prep_w_mat(w)

  # Estimate nuisance parameters
  if (isFALSE(quiet)) print("Propensity score")
  em = nuisance_dss_e(ml_w,wm,x,cf_mat,cv=nfolds,path=path,quiet=quiet)
  if (isFALSE(quiet)) print("Outcome nuisance")
  mm = nuisance_dss_m(ml_y,y,wm,x,cf_mat,cv=nfolds,weights=FALSE,path=path,quiet=quiet)

  # Get the DR scores
  APO = APO_dml(y,mm,wm,em,cf_mat)
  ATE = ATE_dml(APO)

  # Initialize to store results
  cates = NULL

  # Loop over scores and run DR- and NDR-learner
  if (isFALSE(quiet)) print("(N)DR-learner")
  pos = 1
  for (i in 1:(ncol(APO$w_mat)-1)) {
    loc = i+1
    if (isFALSE(compare_all) & i > 1) break
    for (j in loc:(ncol(APO$w_mat))) {
      if (isFALSE(quiet)) print(paste("(N)DR-learner for comparison",toString(i),"and",toString(j)))
      delta = APO$gamma[,j] - APO$gamma[,i]
      if (!is.null(path)) {
        dir.create(paste0(path,"/Comparison",toString(i),toString(j)))
        path_temp = paste0(path,"/Comparison",toString(i),toString(j),"/")
      }
      cates[[pos]] = ndr_core(ml_tau,delta,y,x,APO$w_mat[,c(i,j)],APO$m_mat[,c(i,j)],
                         APO$e_mat[,c(i,j)],cf_mat,xnew=xnew,nfolds=nfolds,path=path_temp,quiet=quiet)
      pos = pos+1
    }
  }
  if (isFALSE(compare_all)) names(cates) = colnames(ATE$delta)[1:length(cates)]
  if (isTRUE(compare_all)) names(cates) = colnames(ATE$delta)
  if (ncol(APO$w_mat) == 1) cates = cates[[1]]

  list("cates"=cates,"APO"=APO,"ATE"=ATE)
}


#' Core function of \code{\link{ndr_learner}}. It executes steps 2 and 4 of algorithms 1 and 2 in Knaus (2020).
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation.
#' Methods can be created by \code{\link{create_method}}.
#' @param delta vector of doubly robust score. E.g. create via \code{\link{ATE_dml}}.
#' @param y Vector of variable to be predicted.
#' @param x Matrix of covariates.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param m_mat n x T+1 matrix with fitted outcome values.
#' @param e_mat n x T+1 matrix with propensity scores.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param xnew Covariate matrix of test sample
#' @param nfolds Number of folds used in cross-validation of ensemble weights (default \code{nfolds=5})
#' @param path Optional path to save the \code{\link{ensemble}} objects for later processing.
#' Saved as Ensemble_Yi where i is the number of the treatment in multiple treatment settings.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return \code{np_grf} returns the grf object(s) and a n x 1 matrix of nuisance parameters
#'
#' @references
#' \itemize{
#' \item Knaus, M. C. (2020). Double machine learning based program evaluation under unconfoundedness.
#'   arXiv preprint arXiv:2003.03191.\url{http://arxiv.org/abs/2003.03191}
#' }
#'
#' @export
#'
ndr_core = function(ml,delta,y,x,w_mat,m_mat,e_mat,cf_mat,
                      xnew=NULL,
                      nfolds=5,
                      path=NULL,
                      quiet=TRUE) {

  if (ncol(cf_mat) != 3) stop("Please provide cf_mat with three columns that was used to create
                              the nuisance parameters and delta.")

  # Predict in sample if no out-of-sample x provided
  if (is.null(xnew)) xnew = x

  # Intialize CATE Matrix
  cates = matrix(0,nrow(xnew),2)
  colnames(cates) = c("DR-Learner","NDR-Learner")

  # Loop over the three folds and calculate DR- and NDR-learner and average
  for (j in 1:3) {
    if (ncol(cf_mat) != 3) stop("Please provide cf_mat with three columns that was used to create the nuisance parameters and delta.")
    ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,j],],y=delta[cf_mat[,j]],xnew=xnew,nfolds=nfolds,weights=TRUE,quiet=quiet)))
    cates[,1] = cates[,1] + 1/3 * ens$ensemble
    cates[,2] = cates[,2] + 1/3 * norm_drl_rcpp(ens$weights,m_mat[cf_mat[,j],],y[cf_mat[,j]],
                                                1*w_mat[cf_mat[,j],],e_mat[cf_mat[,j],])
    if (!is.null(path)) {
      save(ens, file = paste0(path,"Ensemble_DRL_fold",toString(j),".RData"))
    }
    rm(ens); gc()
  }
  return(cates)
}


#' Double sample splitting predictions of propensity score with \code{\link{ensemble}}.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param x Matrix of covariates (n x p matrix)
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param cv Number of cross-validation when estimating ensemble (default 5).
#' @param path Optional path to save the \code{\link{ensemble}} objects for later processing.
#' Saved as Ensemble_Yi where i is the number of the treatment in multiple treatment settings.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return \code{np_grf} returns the grf object(s) and a n x 1 matrix of nuisance parameters
#'
#' @import grf stats
#'
#' @export
#'
nuisance_dss_e = function(ml,w_mat,x,cf_mat,
                       cv=5,
                       path=NULL,
                       quiet=TRUE) {
  # Checks
  if (!is.logical(cf_mat)) stop("Please provide cf_mat as logical matrix. E.g. use function prep_cf_mat")
  if (ncol(cf_mat) != 3) stop("Logical matrix with three columns required for double sample splitting.")
  if (nrow(cf_mat) != nrow(x)) stop("cf_mat matrix nrows different from # of obs.")
  if (nrow(x) != sum(cf_mat)) stop("cf_mat matrix does not sum to number of observations.")

  # Initialize nuisance matrix
  e_mat = matrix(NA,nrow(w_mat),ncol(w_mat))
  colnames(e_mat) = colnames(w_mat)

  ## Calculate treatment probabilities
  # Binary treatment
  if (ncol(w_mat) == 2) {
    # Pscore 1 to 2, 2 to 3, 3 to 1
    ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,1],],y=w_mat[cf_mat[,1],1],
                                  xnew=x[cf_mat[,2],],nfolds=cv,quiet=quiet)))
    e_mat[cf_mat[,2],1] = ens$ensemble
    if (!is.null(path)) {
      save(ens, file = paste0(path,"Ensemble_W_12.RData"))
    }
    rm(ens); gc()
    ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,2],],y=w_mat[cf_mat[,2],1],
                                  xnew=x[cf_mat[,3],],nfolds=cv,quiet=quiet)))
    e_mat[cf_mat[,3],1] = ens$ensemble
    if (!is.null(path)) {
      save(ens, file = paste0(path,"Ensemble_W_23.RData"))
    }
    rm(ens); gc()
    ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,3],],y=w_mat[cf_mat[,3],1],
                                  xnew=x[cf_mat[,1],],nfolds=cv,quiet=quiet)))
    e_mat[cf_mat[,1],1] = ens$ensemble
    if (!is.null(path)) {
      save(ens, file = paste0(path,"Ensemble_W_31.RData"))
    }
    rm(ens); gc()
    e_mat[,2] = 1 - e_mat[,1]
  }
  # Multiple treatment
  else if (ncol(w_mat) > 2) {
    for (i in 1:ncol(w_mat)) {
      # Pscore 1 to 2, 2 to 3, 3 to 1
      ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,1],],y=w_mat[cf_mat[,1],i],
                                    xnew=x[cf_mat[,2],],nfolds=cv,quiet=quiet)))
      e_mat[cf_mat[,2],i] = ens$ensemble
      if (!is.null(path)) {
        save(ens, file = paste0(path,"Ensemble_W",toString(i),"_12.RData"))
      }
      rm(ens); gc()
      ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,2],],y=w_mat[cf_mat[,2],i],
                                    xnew=x[cf_mat[,3],],nfolds=cv,quiet=quiet)))
      e_mat[cf_mat[,3],i] = ens$ensemble
      if (!is.null(path)) {
        save(ens, file = paste0(path,"Ensemble_W",toString(i),"_23.RData"))
      }
      rm(ens); gc()
      ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,3],],y=w_mat[cf_mat[,3],i],
                                    xnew=x[cf_mat[,1],],nfolds=cv,quiet=quiet)))
      e_mat[cf_mat[,1],i] = ens$ensemble
      if (!is.null(path)) {
        save(ens, file = paste0(path,"Ensemble_W",toString(i),"_31.RData"))
      }
      rm(ens); gc()
    }
    e_mat = e_mat / rowSums(e_mat)
  }
  else stop("Provide treatment indicator matrix with at least 2 columns")

  return(e_mat)
}



#' Double sample splitting predictions of outcome nuisance with \code{\link{ensemble}}.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}.
#' @param y Vector of outcome values
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param x Matrix of covariates (n x p matrix)
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param cv Number of cross-validation when estimating ensemble (default 5).
#' @param weights If TRUE, prediction weights of the outcome nuisance extracted and saved (requires to provide a path).
#' @param path Optional path to save the \code{\link{ensemble}} objects for later processing.
#' Saved as Ensemble_Yi where i is the number of the treatment in multiple treatment settings.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return \code{np_grf} returns the grf object(s) and a n x 1 matrix of nuisance parameters
#'
#' @import grf stats
#'
#' @export
#'
nuisance_dss_m = function(ml,y,w_mat,x,cf_mat,
                          cv=5,
                          weights=FALSE,
                          path=NULL,
                          quiet=TRUE) {
  # Checks
  if (!is.logical(cf_mat)) stop("Please provide cf_mat as logical matrix. E.g. use function prep_cf_mat")
  if (ncol(cf_mat) != 3) stop("Logical matrix with three columns required for double sample splitting.")
  if (nrow(cf_mat) != nrow(x)) stop("cf_mat matrix nrows different from # of obs.")
  if (nrow(x) != sum(cf_mat)) stop("cf_mat matrix does not sum to number of observations.")
  if (isTRUE(weights) & is.null(path)) stop("Provide path if weights=TRUE to save ensemble objects with weights for later processing or set weights=FALSE.")

  # Initialize nuisance matrix
  m_mat = matrix(NA,nrow(w_mat),ncol(w_mat))
  colnames(m_mat) = colnames(w_mat)

  ## Calculate outcome predictions
  for (i in 1:ncol(w_mat)) {
    # 1 to 3, 2 to 1, 3 to 2
    ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,1] & w_mat[,i],],y=y[cf_mat[,1] & w_mat[,i]],
                                  xnew=x[cf_mat[,3],],nfolds=cv,weights=weights,quiet=quiet)))
    m_mat[cf_mat[,3],i] = ens$ensemble
    if (!is.null(path)) {
      save(ens, file = paste0(path,"Ensemble_Y",toString(i),"_13.RData"))
    }
    rm(ens); gc()
    ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,2] & w_mat[,i],],y=y[cf_mat[,2] & w_mat[,i]],
                                  xnew=x[cf_mat[,1],],nfolds=cv,weights=weights,quiet=quiet)))
    m_mat[cf_mat[,1],i] = ens$ensemble
    if (!is.null(path)) {
      save(ens, file = paste0(path,"Ensemble_Y",toString(i),"_21.RData"))
    }
    rm(ens); gc()
    ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,3] & w_mat[,i],],y=y[cf_mat[,3] & w_mat[,i]],
                                  xnew=x[cf_mat[,2],],nfolds=cv,weights=weights,quiet=quiet)))
    m_mat[cf_mat[,2],i] = ens$ensemble
    if (!is.null(path)) {
      save(ens, file = paste0(path,"Ensemble_Y",toString(i),"_32.RData"))
    }
    rm(ens); gc()
  }

  return(m_mat)
}


#' This is a wrapper of the C++ function that normalizes the weights for the NDR-learner
#'
#' @param alpha Sparse matrix containing the weights for each observation to be predicted
#' @param m_mat n x T+1 matrix with fitted outcome values.
#' @param y Vector of variable to be predicted.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param e_mat n x T+1 matrix with propensity scores.
#'
#' @return NDR-learner CATEs.
#'
#' @export
#'
norm_drl = function(alpha,m_mat,y,w_mat,e_mat) {
  cates_norm = norm_drl_rcpp(alpha,m_mat,y,1*w_mat,e_mat)
  return(cates_norm)
}

