#' This function estimates the average potential outcomes and average treatment
#' effects using Double Machine Learning (DML).
#'
#' @param y Numeric vector containing the outcome variable.
#' @param w Treatment vector. Provide as factor to control ordering of the treatments,
#' otherwise program orders treatments in ascending order or alphabetically.
#' @param x Covariate matrix.
#' @param ml_w List of methods to be used in ensemble estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param ml_y List of methods to be used in ensemble estimation of outcome regression.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param cf Number of cross-fitting folds for DML (default 5).
#' @param cv Number of cross-validation folds when estimating ensemble if more than one method is defined
#' in \code{ml_w} and/or \code{ml_y} (default 5).
#' @param cl If not NULL, vector with cluster variables
#' @param norm Controls normalization of IPW weights. 0: no normalization, 1: overall normalization,
#' 2: normalization in each cross-fitting fold separately (default).
#' @param weights If TRUE, prediction weights of the outcome nuisance extracted and saved (requires to provide a path).
#' @param path Optional path to save the \code{\link{ensemble}} objects of each cross-fit for later inspection.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#'
#' @return List of an \code{\link{APO_dml}} and an \code{\link{ATE_dml}} object.
#'
#' @export
#'
causalDML = function(y,w,x,
                     ml_w = list(create_method("forest_grf")),
                     ml_y = list(create_method("forest_grf")),
                     cf=5,
                     cv=5,
                     cl=NULL,
                     norm=2,
                     weights=FALSE,
                     path = NULL,
                     quiet=TRUE) {

  # Create important matrices
  cfm = prep_cf_mat(length(y),cf)
  wm = prep_w_mat(w)

  # Estimate pscore
  em = nuisance_e(ml_w,wm,x,cfm,cv=cv,path=path,quiet=quiet)

  # Estimate outcomes
  mm = nuisance_m(ml_y,y,wm,x,cfm,cv=cv,weights=weights,path=path,quiet=quiet)

  # APOs
  APO = APO_dml(y,mm,wm,em,norm=norm,cf_mat=cfm)

  # ATE
  ATE = ATE_dml(APO)

  list("APO"=APO,"ATE"=ATE)
}


#' This function estimates the average potential outcomes and average treatment
#' effects using Double Machine Learning (DML).
#'
#' More recent version of \code{\link{causalDML}} with more functions and especially
#' a more precise function name.
#'
#' @param y Numeric vector containing the outcome variable.
#' @param w Treatment vector. Provide as factor to control ordering of the treatments,
#' otherwise program orders treatments in ascending order or alphabetically.
#' @param x Covariate matrix.
#' @param ml_w List of methods to be used in ensemble estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param ml_y List of methods to be used in ensemble estimation of outcome regression.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param cf Number of cross-fitting folds for DML (default 5).
#' @param cv Number of cross-validation folds when estimating ensemble if more than one method is defined
#' in \code{ml_w} and/or \code{ml_y} (default 5).
#' @param cl If not NULL, vector with cluster variables
#' @param norm Controls normalization of IPW weights. 0: no normalization, 1: overall normalization,
#' 2: normalization in each cross-fitting fold separately (default).
#' @param weights If TRUE, prediction weights of the outcome nuisance extracted and saved (requires to provide a path).
#' @param path Optional path to save the \code{\link{ensemble}} objects of each cross-fit for later inspection.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#' @param e_mat Optional n x T+1 matrix with propensity scores calculated outside of function.
#' @param m_mat Optional n x T+1 matrix fitted outcome values calculated outside of function.
#' @param cf_mat Optional prespecified logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#'
#' @return List of an \code{\link{APO_dml}} and an \code{\link{ATE_dml}} object.
#'
#' @export
#'
DML_aipw = function(y,w,x,
                     ml_w = list(create_method("forest_grf")),
                     ml_y = list(create_method("forest_grf")),
                     cf=5,
                     cv=5,
                     cl=NULL,
                     norm=2,
                     weights=FALSE,
                     path = NULL,
                     quiet=TRUE,
                     e_mat=NULL,
                     m_mat=NULL,
                     cf_mat=NULL) {
  # Create important matrices
  wm = prep_w_mat(w)
  if (is.null(cf_mat)) cfm = prep_cf_mat(length(y),cf,wm,cl)
  else cfm = cf_mat

  # Estimate pscore
  if(is.null(e_mat)) em = nuisance_e(ml_w,wm,x,cfm,cv=cv,path=path,quiet=quiet)
  else em = e_mat

  # Estimate outcomes
  if(is.null(m_mat)) mm = nuisance_m(ml_y,y,wm,x,cfm,cv=cv,weights=weights,path=path,quiet=quiet)
  else mm = m_mat

  # APOs
  APO = APO_dml(y,mm,wm,em,norm=norm,cf_mat=cfm)

  # ATE
  ATE = ATE_dml(APO)

  list("APO"=APO,"ATE"=ATE)
}


#' This function estimates the average potential outcomes for all treatment levels.
#'
#' @param y A numeric vector containing the outcome variable.
#' @param m_mat n x T+1 matrix with fitted outcome values.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param e_mat n x T+1 matrix with propensity scores.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param norm Controls normalization of IPW weights. \code{norm = 1)}: overall normalization, \code{norm = 2)}: normalization in each
#'  cross-fitting fold separately (default), otherwise no normalization.
#' @param cl If not NULL, vector with cluster variables.
#'
#' @return Returns an \code{APO_dml} object:
#'          \item{results}{Point estimate and standard error of APOs.}
#'          \item{gamma}{n x T+1 matrix containing the APO DR scores.}
#'          \item{y}{Vector of outcomes used in estimation.}
#'          \item{m_mat}{Matrix of outcome predictions used in estimation.}
#'          \item{w_mat}{Matrix of treatment indicators (n x T+1) used in estimation.}
#'          \item{e_mat}{Matrix of propensity score used in estimation.}
#'          \item{cf_mat}{Matrix with k columns of indicators representing the different folds used in estimation.}
#'          \item{cl}{Cluster variable if specified.}
#'
#' @import matrixStats
#'
#' @export
#'
APO_dml = function(y,
                   m_mat,
                   w_mat,
                   e_mat,
                   cf_mat,
                   norm=2,
                   cl=NULL) {
  # Initialize matrices
  lambda_ipw = matrix(0,nrow(w_mat),ncol(w_mat))
  gamma_mat = matrix(NA,nrow(w_mat),ncol(w_mat))
  res = matrix(NA,ncol(w_mat),2)
  rownames(res) = colnames(w_mat)
  colnames(res) = c("APO","SE")

  # Loop over all treatments to get p-score
  for (i in 1:ncol(w_mat)) {
    lambda_ipw[,i] = w_mat[,i] / e_mat[,i]
    if (norm == 1)  lambda_ipw[,i] = lambda_ipw[,i] / sum(lambda_ipw[,i]) * nrow(w_mat)
    if (norm == 2) {
      for (j in 1:ncol(cf_mat)) {
        lambda_ipw[cf_mat[,j],i] = lambda_ipw[cf_mat[,j],i] /
                                        sum(lambda_ipw[cf_mat[,j],i]) * nrow(w_mat[cf_mat[,j],])
      }
    }
  }

  # Calculate APO DR score
  for (i in 1:ncol(w_mat)) {
    gamma_mat[,i] = m_mat[,i] + lambda_ipw[,i] * (y - m_mat[,i])
  }

  # Calculate APO
  res[,1] = colMeans(gamma_mat)

  # Calculate SE for APO
  if (is.null(cl)) {
    for (i in 1:ncol(w_mat)) {
      res[i,2] = sqrt(mean((gamma_mat[,i] - mean(gamma_mat[,i]))^2) / nrow(w_mat))
    }
  } else {
    for (i in 1:ncol(w_mat)) {
      res[i,2] = sqrt(sum(tapply(gamma_mat[,i] - mean(gamma_mat[,i]), cl, sum)^2) / nrow(w_mat)^2)
    }
  }

  output = list("results"=res,"gamma"=gamma_mat,"y"=y,"m_mat"=m_mat,"w_mat"=w_mat,"e_mat"=e_mat,"cf_mat"=cf_mat,"cl"=cl)
  class(output) = "APO_dml"
  output
}


#' This function calculates the average treatment effects (on the treated) for all combinations
#' of potentially multiple treatments.
#'
#' @param APO Object of class \code{\link{APO_dml}} or \code{\link{APO_dml_atet}}.
#'
#' @return Returns an \code{ATE_dml} object:
#'          \item{results}{Point estimate and standard error of APOs.}
#'          \item{delta}{n x T+1 matrix containing the ATE DR scores.}
#'
#' @export
#'
ATE_dml = function(APO) {

  # Initialize matrices
  delta = matrix(NA,nrow(APO$w_mat),sum(1:(ncol(APO$w_mat)-1)))
  res = matrix(NA,sum(1:(ncol(APO$w_mat)-1)),4)
  if (class(APO) == "APO_dml") colnames(res) = c("ATE","SE","t","p")
  if (class(APO) == "APO_dml_atet") colnames(res) = c("ATET","SE","t","p")
  rownames(res) = rep("Platzhalter",nrow(res))

  # Calcuate ATE for all combinations
  pos = 1
  for (i in 1:(ncol(APO$w_mat)-1)) {
    loc = i+1
    for (j in loc:(ncol(APO$w_mat))) {
      delta[,pos] = APO$gamma[,j] - APO$gamma[,i]
      res[pos,1] = mean(delta[,pos])
      if (is.null(APO$cl)) {
        res[pos,2] = sqrt(mean((delta[,pos]-mean(delta[,pos]))^2) / nrow(APO$w_mat))
      } else {
        res[pos,2] =  sqrt(sum(tapply(delta[,pos] - mean(delta[,pos]), APO$cl, sum)^2) / nrow(APO$w_mat)^2)
      }
      rownames(res)[pos] = paste(colnames(APO$w_mat)[j],"-",colnames(APO$w_mat)[i])
      pos = pos + 1
    }
  }
  colnames(delta) = rownames(res)

  # t-stat
  res[,3] = res[,1] / res[,2]
  # p-value
  res[,4] = 2 * stats::pt(abs(res[,3]),nrow(APO$w_mat),lower = FALSE)

  output = list("results"=res,"delta"=delta)
  class(output) = "ATE_dml"
  output
}



#' This function calculates the average potential outcomes for the average treatment effect on the treated.
#'
#' @param y A numeric vector containing the outcome variable.
#' @param m_mat n x T+1 matrix with fitted outcome values.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param e_mat n x T+1 matrix with (generalized) propensity scores.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param treated An integer describing the column number of the treated in w_mat for which ATET
#' should be calculated (default 2).
#' @param norm Controls normalization of IPW weights. \code{norm = 1)}: overall normalization, \code{norm = 2)}:
#' normalization in each cross-fitting fold separately (default), otherwise no normalization.
#' @param cl If not NULL, vector with cluster variables.
#'
#' @return \code{APO_dmlmt} returns a list containing the \code{results} and a matrix
#'          with the individual values of the doubly robust score.
#' @export
#'
APO_dml_atet = function(y,
                        m_mat,
                        w_mat,
                        e_mat,
                        cf_mat,
                        treated=2,
                        norm=2,
                        cl=NULL) {

  # Initialize matrices
  lambda_ipw = matrix(0,nrow(w_mat),ncol(w_mat))
  gamma_mat = matrix(NA,nrow(w_mat),ncol(w_mat))
  res = matrix(NA,ncol(w_mat),2)
  rownames(res) = colnames(w_mat)
  colnames(res) = c("ATET APO","SE")

  # Get pscores of treated
  ewx = e_mat[,treated]
  ew = mean(w_mat[,treated])

  # Get IPW weights
  for (i in 1:ncol(w_mat)) {
    lambda_ipw[,i] =  ((w_mat[,i]) * ewx) / (ew * e_mat[,i])
    if (norm == 1)  lambda_ipw[,i] = lambda_ipw[,i] / sum(lambda_ipw[,i]) * nrow(w_mat)
    if (norm == 2) {
      for (j in 1:ncol(cf_mat)) {
        lambda_ipw[cf_mat[,j],i] = lambda_ipw[cf_mat[,j],i] /
          sum(lambda_ipw[cf_mat[,j],i]) * nrow(w_mat[cf_mat[,j],])
      }
    }
  }

  # APO DR score for treated
  for (i in 1:ncol(w_mat)) {
    gamma_mat[,i] = ((w_mat[,treated]) * m_mat[,i]) / ew + lambda_ipw[,i] * (y - m_mat[,i])
  }

  # Calculate APO
  res[,1] = colMeans(gamma_mat)

  # Calculate SE for APO
  if (is.null(cl)) {
    for (i in 1:ncol(w_mat)) {
      res[i,2] = sqrt(mean((gamma_mat[,i] - mean(gamma_mat[,i]))^2) / nrow(w_mat))
    }
  } else {
    for (i in 1:ncol(w_mat)) {
      res[i,2] = sqrt(sum(tapply(gamma_mat[,i] - mean(gamma_mat[,i]), cl, sum)^2) / nrow(w_mat)^2)
    }
  }

  output = list("results"=res,"gamma"=gamma_mat,"y"=y,"m_mat"=m_mat,"w_mat"=w_mat,"e_mat"=e_mat,"cl"=cl)
  class(output) = "APO_dml_atet"
  output
}



#' \code{summary} method for class \code{\link{APO_dml}}
#'
#' @param APO_dml Object of class \code{\link{APO_dml}}.
#'
#' @export
#'
summary.APO_dml = function(APO_dml) {
  print(APO_dml$results)
}

#' \code{summary} method for class \code{\link{ATE_dml}}
#'
#' @param ATE_dml Object of class \code{\link{ATE_dml}}.
#' @param ... further argumenst passed to \code{printCoefmat}
#'
#' @export
#'
summary.ATE_dml = function(ATE_dml,...) {
  printCoefmat(ATE_dml$results,has.Pvalue = TRUE)
}


#' \code{plot} method for class \code{\link{APO_dml}}
#'
#' @param APO_dml Object of class \code{\link{APO_dml}}.
#' @param label Optional vector of treatment labels.
#' @param sl Significance level for confidence intervals (default 0.05).
#'
#' @import ggplot2
#'
#' @export
#'
plot.APO_dml = function(APO_dml,label = NULL,sl=0.05) {
  PO = w = cil = ciu = NULL
  qt = qt(1-sl / 2,nrow(APO_dml$gamma)-1)

  if(is.null(label)) label = factor(colnames(APO_dml$w_mat),levels = colnames(APO_dml$w_mat))

  df = data.frame(PO=APO_dml$results[,1],w=factor(label,levels = label),
                  cil=APO_dml$results[,1]-qt*APO_dml$results[,2],
                  ciu=APO_dml$results[,1]+qt*APO_dml$results[,2])

  g = ggplot(df,aes(x=w,y=PO,ymin=cil,ymax=ciu)) + geom_point(size=2.5) +  geom_errorbar(width=0.15)  +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Average potential outcome") + theme(text=element_text(family="serif",size = 16, colour="black"),
                                              axis.title.x=element_blank())
  return(g)
}



#' Double Machine Learning estimation of partially linear model
#'
#' This function estimates the parameter of interest in a partially linear model
#' using the residual-on-residual representation of Robinson (1988)
#' and flexibly estimated nuisance parameters following Chernozhukov et al. (2018).
#'
#' @param y Numeric vector containing the outcome variable.
#' @param w Treatment vector (binary or continuous).
#' @param x Covariate matrix.
#' @param ml_w List of methods to be used in ensemble estimation of treatment regression.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param ml_y List of methods to be used in ensemble estimation of outcome regression.
#' Methods can be created by \code{\link{create_method}}. Default is an untuned honest
#' \code{\link{regression_forest}}.
#' @param cf Number of cross-fitting folds for DML (default 5).
#' @param cv Number of cross-validation folds when estimating ensemble if more than one method is defined
#' in \code{ml_w} and/or \code{ml_y} (default 5).
#' @param weights If TRUE, prediction weights of the outcome nuisance extracted and saved (requires to provide a path).
#' @param path Optional path to save the \code{\link{ensemble}} objects of each cross-fit for later inspection.
#' @param quiet If FALSE, ensemble estimators print method that is currently running.
#' @param e_hat Optional vector of predicted treatment outside of the function.
#' @param m_hat Optional vector of predicted outcome outside of the function.
#' @param cf_mat Optional prespecified logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#'
#' @return Returns a \code{DML_partial_linear} object:
#'          \item{results}{Point estimate, standard error, t- and p-value of estimated effect.}
#'          \item{e_hat}{Predicted treatment}
#'          \item{m_hat}{Predicted outcomes}
#'          \item{w}{Treatment vector used in estimation.}
#'          \item{y}{Vector of outcomes used in estimation.}
#'          \item{cf_mat}{Matrix with k columns of indicators representing the different folds used in estimation.}
#'          \item{path}{Path where results are stored if specified, otherwise NULL.}
#'
#' @references
#' \itemize{
#' \item Robinson, P. M. (1988). Root-N-consistent semiparametric regression. Econometrica, 931-954.
#' \item Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018).
#' Double/Debiased machine learning for treatment and structuralparameters.The Econometrics Journal,21(1), C1-C68
#' }
#'
#' @export
#'
DML_partial_linear = function(y,w,x,
                              ml_w = list(create_method("forest_grf")),
                              ml_y = list(create_method("forest_grf")),
                              cf=5,
                              cv=5,
                              weights=FALSE,
                              path=NULL,
                              quiet=TRUE,
                              e_hat=NULL,
                              m_hat=NULL,
                              cf_mat=NULL) {
  # Get important matrices
  if (is.null(cf_mat)) cfm = prep_cf_mat(length(y),cf)
  else cfm = cf_mat

  # Manage paths
  path_w = path_y = NULL
  if (!is.null(path)) {
    path_w = paste0(path,"W_")
    path_y = paste0(path,"Y_")
  }

  # Estimate treatment model
  n = length(y)
  if(is.null(e_hat)) ex = nuisance_m(ml_w,w,matrix(T,nrow=n,ncol=1),x,cfm,cv=cv,
                                     path=path_w,quiet=quiet)
  else ex = e_hat
  v = w - ex

  # Estimate outcomes
  if(is.null(m_hat)) mx = nuisance_m(ml_y,y,matrix(T,nrow=n,ncol=1),x,cfm,cv=cv,weights=weights,
                                     path=path_y,quiet=quiet)
  else mx = m_hat
  u = y - mx

  # Results
  res = matrix(NA,1,4)
  colnames(res) = c("Coefficient","SE","t","p")
  pl = lm(u ~ 0 + v)
  res[,1] = coef(pl)
  res[,2] = sqrt((1/sum(v^2))^2 * sum((v*pl$residuals)^2))
  res[,3] = res[,1] / res[,2]
  res[,4] = 2 * stats::pt(abs(res[,3]),n,lower = FALSE)

  output = list("result"=res,"e_hat"=ex,"m_hat"=mx,"w"=w,"y"=w,"cf_mat"=cfm,"path"=path)
  class(output) = "DML_partial_linear"
  output
}


#' \code{summary} method for class \code{\link{DML_partial_linear}}
#'
#' @param DML_partial_linear Object of class \code{\link{DML_partial_linear}}.
#' @param ... further arguments passed to \code{printCoefmat}
#'
#' @export
#'
summary.DML_partial_linear = function(DML_partial_linear,...) {
  printCoefmat(DML_partial_linear$result,has.Pvalue = TRUE)
}
