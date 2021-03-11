#' DR-learner
#'
#' This function produces out-of-sample predictions of conditional average treatment effects (CATEs) for
#' all individuals in the sample using the DR-learner (Kennedy, 2020).
#' The involved predictions are based on \code{\link{ensemble}}. See algorithm 1 Knaus (2020) for details.
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
#' for later processing. Outputs of each fold i saved in new subfolder DR_foldi.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return  \item{cates}{Returns a n x number of comparisons matrix (one column in case of binary treatment)
#' containing the predicted CATEs.}
#'          \item{list}{A list of the four \code{\link{dr_oos}} outputs.}
#'          \item{cf_mat}{Matrix with k columns of indicators representing the different folds used in estimation.}
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
dr_learner = function(y,w,x,
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
  cfm = prep_cf_mat(length(y),4,wm)

  # How many comparisons are specified?
  if (isTRUE(compare_all)) num_comp = ncol(wm)*(ncol(wm)-1)/2
  else num_comp = ncol(wm)-1

  # Initialize vectors to store results
  list_dr_oos = vector("list",4)
  cates = matrix(NA,nrow(x),num_comp)
  path_orig = path

  # Loop over four folds
  for (i in 1:4) {
    if (isFALSE(quiet)) print(paste("DR-learner fold", toString(i)))
    oos = cfm[,i]
    if (!is.null(path)) {
      dir.create(paste0(path_orig,"/DR_fold",toString(i)))
      path = paste0(path_orig,"/DR_fold",toString(i),"/")
    }
    list_dr_oos[[i]] = dr_oos(y[!oos],w[!oos],x[!oos,],ml_w = ml_w,
                                ml_y = ml_y,ml_tau = ml_tau,
                                cf_mat = cfm[!oos,-i],compare_all = compare_all,xnew=x[oos,],
                                path=path,quiet=quiet)
    cates[oos,] = list_dr_oos[[i]]$cates
  }

  colnames(cates) = colnames(list_dr_oos[[i]]$cates)

  list("cates"=cates,"list"=list_dr_oos,"cf_mat"=cfm)
}


#' Out-of-sample prediction with DR-learner
#'
#' This function produces out-of-sample predictions of conditional average treatment effects (CATEs) using the
#' DR-learner (Kennedy, 2020).
#' It executes steps 1 to 4 of algorithm 1 in Knaus (2020) as part of the \code{\link{dr_learner}}.
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
#' @return  \item{cates}{Returns a n x number of comparisons matrix (one column in case of binary treatment)
#' containing the predicted CATEs.}
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
dr_oos = function(y,w,x,xnew,
                   ml_w = list(create_method("forest_grf")),
                   ml_y = list(create_method("forest_grf")),
                   ml_tau = list(create_method("forest_grf")),
                   cf_mat = NULL,
                   compare_all = TRUE,
                   nfolds=5,
                   path=NULL,
                   quiet=TRUE) {
  # Create indicator matrices
  wm = prep_w_mat(w)
  if (is.null(cf_mat)) cf_mat = prep_cf_mat(length(y),3,wm)

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
  path_orig = path

  # Loop over scores and run DR-learner
  if (isFALSE(quiet)) print("DR-learner")
  pos = 1
  for (i in 1:(ncol(APO$w_mat)-1)) {
    loc = i+1
    if (isFALSE(compare_all) & i > 1) break
    for (j in loc:(ncol(APO$w_mat))) {
      if (isFALSE(quiet)) print(paste("DR-learner for comparison",toString(j),"and",toString(i)))
      delta = APO$gamma[,j] - APO$gamma[,i]
      if (!is.null(path)) {
        dir.create(paste0(path_orig,"/Comparison",toString(i),toString(j)))
        path_temp = paste0(path_orig,"/Comparison",toString(i),toString(j),"/")
      }
      cates = cbind( cates , dr_core(ml_tau,delta,y,x,APO$w_mat[,c(i,j)],APO$m_mat[,c(i,j)],
                              APO$e_mat[,c(i,j)],cf_mat,xnew=xnew,nfolds=nfolds,path=path,quiet=quiet) )
      pos = pos+1
    }
  }
  if (isFALSE(compare_all)) colnames(cates) = colnames(ATE$delta)[1:ncol(cates)]
  if (isTRUE(compare_all)) colnames(cates) = colnames(ATE$delta)

  list("cates"=cates,"APO"=APO,"ATE"=ATE)
}


#' Core function of \code{\link{dr_learner}}. It executes steps 2 and 4 of algorithm 1 in Knaus (2020).
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
#' @return Vector with predicted CATEs.
#'
#' @references
#' \itemize{
#' \item Knaus, M. C. (2020). Double machine learning based program evaluation under unconfoundedness.
#'   arXiv preprint arXiv:2003.03191.\url{http://arxiv.org/abs/2003.03191}
#' }
#'
#' @export
#'
dr_core = function(ml,delta,y,x,w_mat,m_mat,e_mat,cf_mat,
                    xnew=NULL,
                    nfolds=5,
                    path=NULL,
                    quiet=TRUE) {

  if (ncol(cf_mat) != 3) stop("Please provide cf_mat with three columns that was used to create
                              the nuisance parameters and delta.")

  # Predict in sample if no out-of-sample x provided
  if (is.null(xnew)) xnew = x

  # Intialize CATE Matrix
  cate = rep(0,nrow(xnew))

  # Loop over the three folds and calculate DR-learner and average
  for (j in 1:3) {
    if (ncol(cf_mat) != 3) stop("Please provide cf_mat with three columns that was used to create the nuisance parameters and delta.")
    ens = do.call(ensemble,c(list(ml=ml,x=x[cf_mat[,j],],y=delta[cf_mat[,j]],xnew=xnew,nfolds=nfolds,quiet=quiet)))
    cate = cate + 1/3 * ens$ensemble
    if (!is.null(path)) {
      save(ens, file = paste0(path,"Ensemble_DRL_fold",toString(j),".RData"))
    }
    rm(ens); gc()
  }
  return(cate)
}

