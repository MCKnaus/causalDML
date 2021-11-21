#' Create interactions and polynomials
#'
#' This function takes a dataset and strings with variable names as input and creates interactions and polynomials
#'
#' @param data Matrix with the main variables.
#' @param int Vector of strings with the variables to be interacted.
#' @param int_d Degree of interactions created. Default is 2.
#' @param poly Vector of strings with the variables for which polynomials should be created.
#' @param poly_d Degree of polynomials to be created. Default is 2.
#' @param log Vector of strings with the variables for which the logged versions should be added.
#'
#' @return Matrix including the main variables and the newly generated ones.#'
#'
#' @export
#'
design_matrix = function(data,int=NULL,int_d=2,poly=NULL,poly_d=2,log=NULL) {

  # Part for interactions
  if (!is.null(int)) {
    int = paste0("+(",paste0(int,collapse="+"),")^",toString(int_d))
  }

  # Part for polynomials
  if (!is.null(poly)) {
    poly = paste0("+poly(",paste0(poly,collapse=paste0(",",toString(poly_d),",raw=TRUE)+poly(")),",",toString(poly_d),",raw=TRUE)",
                   "-(",paste0(paste0(poly,collapse="+")),")")
  }

  # Part for logs
  if (!is.null(log)) {
    # Check whether some variables can't be logged because not positive
    ind_neg = colMins(data[,log])<=0
    if (sum(ind_neg)>0) {
      cat("\n Following variables not modified to be logged because of non-positive values:",paste(colnames(data[,log])[ind_neg]),"\n" )
      log = log[!ind_neg]
    }
    log = paste0("+log(",paste0(log,collapse=")+log("),")")
  }

  # Combine the three parts
  fmla = as.formula(paste("~0",int,poly,log))

  # Generate matrix
  data = model.matrix(fmla,data=as.data.frame(data))
  # Clean variable names to make sense
  colnames(data) = gsub("poly\\(","",colnames(data))
  colnames(data) = gsub(paste0(", ",toString(poly_d),", raw = TRUE)"),"",colnames(data))
  colnames(data) = gsub("log\\(","ln_",colnames(data))
  colnames(data) = gsub("\\)","",colnames(data))

  return(data)
}


#' Data screening
#'
#' This function takes a matrix of data and removes
#' 1. Variables without variation
#' 2. Dummy variables where one group is nearly empty (optional in one of both treatment groups)
#' 3. Redundant (highly correlated variables)
#'
#' @param data Matrix the variables to be screened.
#' @param treat Optional binary treatment vector if screening should be done within treatment groups
#' @param bin_cut Cut-off fraction under which nearly empty binary variables should be removed. Default 0.01.
#' @param corr_cut Cut-off above which highly correlated varialbes should be removed. Default 0.99.
#' @param print Shows details about the reomved variables at each step
#'
#' @return Screened matrix
#'
#' @import matrixStats
#'
#' @export
#'
data_screen = function(data,treat=NULL,bin_cut=0.01,corr_cut=0.99,print=FALSE) {

  ## Kick out variables with no variation
  # Identify the names
  nm_del = colnames(data)[colSds(data) == 0]
  # Describe identified variables
  if (print==TRUE) {
    cat("\n\n Variables with no variation:",nm_del,"\n\n")
  }
  # Remove identified variables
  if (identical(nm_del, character(0)) == FALSE) data = data[,!colnames(data) %in% nm_del]

  ## Remove dummy variables lower than threshold in one of the two treatment groups
  # Identify dummies
  bin = apply(data,2,function(x) { all(x %in% 0:1) })

  # Calculate means of all variables and check whether they are potentially close to 0 or 1
  if (is.null(treat)) {
    mean = colMeans(data)
    bel_cut = (mean<bin_cut | mean > (1-bin_cut))
  } else {
    mean1 = colMeans(data[d==1,])
    mean0 = colMeans(data[d==0,])
    bel_cut = (mean1<bin_cut | mean1 > (1-bin_cut) | mean0<bin_cut | mean0 > (1-bin_cut))
  }

  # Identify names that are binary and close to 0 and 1
  nm_del = colnames(data)[bin & bel_cut]
  if (print==TRUE) {
    cat("\n\n Dummy variables close to 0 or 1:",nm_del,"\n\n")
  }

  # Remove identified variables
  if (identical(nm_del, character(0)) == FALSE) data = data[,!colnames(data) %in% nm_del]

  ## Remove all redundant (nearly perfectly correlated) variables
  # Calculate correlation matrix and consider only upper diagonal
  cor = (abs(cor(data))>corr_cut)
  cor[lower.tri(cor, diag=TRUE)] = FALSE

  # Identify names of redundant variables
  nm_del = colnames(cor)[colSums(cor)>0]

  if (print==TRUE) {
    cat("\n\n Variables (nearly) perfectly correlated:",nm_del,"\n\n")
  }
  if (identical(nm_del, character(0)) == FALSE) data = data[,colSums(cor)==0]

  return(data)
}
