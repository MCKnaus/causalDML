#' Implementation of Heiler & Knaus (2025) decomposition. 
#'
#' @param Y Numeric vector containing the outcome variable.
#' @param A Aggregate treatment vector, e.g. binary. Program finds mapping to effective
#' treatment automatically if possible.
#' @param G Heterogeneity group vector. Provide as factor to control ordering,
#' otherwise program orders treatments in ascending order or alphabetically.
#' @param T_mat Logical matrix of effective treatment indicators (n x J). 
#' For example created by \code{\link{prep_w_mat}}.
#' @param e_mat n x J matrix with propensity scores.
#' @param m_mat n x J matrix with fitted outcome values.
#' @param sampling_weights Optional vector of sampling weights.
#' @param cl If not NULL, vector with cluster variables.
#'
#' @return Returns an \code{HK2_decomposition} object:
#'          \item{parameter}{14 x # of heterogeneity groups x # of treatment aggregates x 2 array
#'          storing point estimates and standard error of target and intermediate parameters. The ordering is
#'          c("GM","AGM","d0","s1","s2","s3","d1","d2","d3","Cov(etX,mut|Xg)","d4","SRCT2","d4'","d5").}
#'          \item{IFs}{14 x # of heterogeneity groups x # of treatment aggregates x n x 6 array
#'          storing the influence fcts and its components corresponding to each parameter for further use.}
#'          \item{mapping}{Logical matrix storing the mapping of aggregate and effective treatment.}
#'          \item{label}{List of labels for further useage.}
#'          \item{cl}{Cluster variable if specified.}
#'
#' @references
#' \itemize{
#' \item Heiler, P., Knaus, M.C. (2025). Heterogeneity analysis with heterogeneous treatments.
#' }
#'
#' @export
#'
HK2_decomposition = function(Y,A,G,
                            T_mat,
                            e_mat,
                            m_mat,
                            sampling_weights = NULL,
                            cl = NULL) {
  # Sanity check
  if (!(identical(colnames(e_mat), colnames(m_mat)) &&
        identical(colnames(m_mat), colnames(T_mat)))) {
    warning("Column names of e_mat, m_mat, and T_mat are not identical. Make sure this is not an issue.")
  }
  
  ### Define crucial ingredients
  n = length(Y)
  J = ncol(T_mat)
  num_A = length(table(A)) + 1
  num_G = length(table(G)) + 1
  ones = rep(1,n)
  A_mat = prep_w_mat(A)
  G_mat = prep_w_mat(unlist(G))
  
  A_label = c("All",colnames(A_mat))
  T_label = colnames(T_mat)
  G_label = c("Unconditional",colnames(G_mat))
  
  treatment_mapping = cbind(rep(TRUE,J), # all in one aggregation
                            treatment_aggregation_mapping(A_mat,T_mat)) #recover the implied agg
  colnames(treatment_mapping) = A_label
  
  #########################################################
  ### Prepare i-varying (raw) nuisance parameters
  # ti varying
  NP_ti_varying = array(NA, dim = c(3,n,J))
  dimnames(NP_ti_varying)[[1]] = c("1(T = t)","et(x)","μt(x)")
  dimnames(NP_ti_varying)[[3]] = T_label
  NP_ti_varying[1,,] = T_mat
  NP_ti_varying[2,,] = e_mat
  NP_ti_varying[3,,] = m_mat
  
  # ai varying (part 1)
  NP_ai_varying = array(NA, dim = c(3,n,num_A))
  dimnames(NP_ai_varying)[[1]] = c("1(T ∈ Ta)","ea(x)","μa(x)")
  dimnames(NP_ai_varying)[[3]] = A_label
  NP_ai_varying[1,,] = cbind(ones,A_mat)
  NP_ai_varying[2,,] = NP_ti_varying[2,,] %*% treatment_mapping
  
  # tai varying
  NP_tai_varying = array(0, dim = c(1,num_A,n,J))
  dimnames(NP_tai_varying)[[1]] = "eta(x)"
  dimnames(NP_tai_varying)[[2]] = A_label
  dimnames(NP_tai_varying)[[4]] = T_label
  for (a in 1:num_A) {
    NP_tai_varying[1,a,,treatment_mapping[,a]] = NP_ti_varying[2,,treatment_mapping[,a]] / NP_ai_varying[2,,a]
  }
  
  # ai varying (part 2)
  for (a in 1:num_A) {
    NP_ai_varying[3,,a] = rowSums(NP_tai_varying[1,a,,] * NP_ti_varying[3,,])
  }
  
  # gi varying
  NP_gi_varying = array(NA, dim = c(1,n,num_G))
  dimnames(NP_gi_varying)[[1]] = "1(X ∈ Xg)"
  dimnames(NP_gi_varying)[[3]] = G_label
  NP_gi_varying[1,,] = cbind(ones,G_mat)
  
  # agi varying
  NP_agi_varying = array(NA, dim = c(1,num_A,n,num_G))
  dimnames(NP_agi_varying)[[1]] = "1(T ∈ Ta,X ∈ Xg)"
  dimnames(NP_agi_varying)[[2]] = A_label
  dimnames(NP_agi_varying)[[4]] = G_label
  for (a in 1:num_A) {
    NP_agi_varying[1,a,,] = NP_gi_varying[1,,] * NP_ai_varying[1,,a]  
  }
  
  #########################################################
  ### Prepare aggregated nuisance parameters and their IFs
  # a varying
  NP_a_varying = array(NA, dim = c(1,num_A,2))
  NP_IF_a_varying = array(NA, dim = c(1,num_A,n,6))
  dimnames(NP_a_varying)[[1]] = dimnames(NP_IF_a_varying)[[1]] = "ea"
  dimnames(NP_a_varying)[[2]] = dimnames(NP_IF_a_varying)[[2]] = A_label
  for (a in 1:num_A) {
    # Unconditional prob e_a
    temp = IF_maker(ones,NP_ai_varying[1,,a],-ones,
                    sampling_weights,cl)
    NP_a_varying[1,a,1] = temp$estimate
    NP_a_varying[1,a,2] = temp$se
    NP_IF_a_varying[1,a,,] = temp$IF
  }
  
  # g varying
  NP_g_varying = array(NA, dim = c(1,num_G,2))
  NP_IF_g_varying = array(NA, dim = c(1,num_G,n,6))
  dimnames(NP_g_varying)[[1]] = dimnames(NP_IF_g_varying)[[1]] = "eg"
  dimnames(NP_g_varying)[[2]] = dimnames(NP_IF_g_varying)[[2]] = G_label
  for (g in 1:num_G) {
    # Unconditional prob e_g
    temp = IF_maker(ones,NP_gi_varying[1,,g],-ones,
                    sampling_weights,cl)
    NP_g_varying[1,g,1] = temp$estimate
    NP_g_varying[1,g,2] = temp$se
    NP_IF_g_varying[1,g,,] = temp$IF
  }
  
  # t varying
  NP_t_varying = array(NA, dim = c(2,J,2))
  NP_IF_t_varying = array(NA, dim = c(2,J,n,6))
  dimnames(NP_t_varying)[[1]] = dimnames(NP_IF_t_varying)[[1]] = c("et","μt")
  dimnames(NP_t_varying)[[2]] = dimnames(NP_IF_t_varying)[[2]] = T_label
  for (t in 1:J) {
    # e_t
    temp = IF_maker(ones,NP_ti_varying[1,,t],-ones,
                    sampling_weights,cl)
    NP_t_varying[1,t,1] = temp$estimate
    NP_t_varying[1,t,2] = temp$se
    NP_IF_t_varying[1,t,,] = temp$IF
    # mu_t
    temp = IF_maker(ones,
                    NP_ti_varying[3,,t] #mhat
                    + NP_ti_varying[1,,t] * (Y - NP_ti_varying[3,,t]) / NP_ti_varying[2,,t],
                    -ones,
                    sampling_weights,cl)
    NP_t_varying[2,t,1] = temp$estimate
    NP_t_varying[2,t,2] = temp$se
    NP_IF_t_varying[2,t,,] = temp$IF
  }
  
  # ta varying
  NP_ta_varying = array(NA, dim = c(1,num_A,J,2))
  NP_IF_ta_varying = array(NA, dim = c(1,num_A,J,n,6))
  dimnames(NP_ta_varying)[[1]] = dimnames(NP_IF_ta_varying)[[1]] = "eta"
  dimnames(NP_ta_varying)[[2]] = dimnames(NP_IF_ta_varying)[[2]] = A_label
  dimnames(NP_ta_varying)[[3]] = dimnames(NP_IF_ta_varying)[[3]] = T_label
  for (a in 1:num_A) {
    for (t in 1:J) {
      temp = IF_maker(NP_ai_varying[1,,a] / NP_a_varying[1,a,1],
                      NP_ti_varying[1,,t],
                      -ones,
                      sampling_weights,cl)
      NP_ta_varying[1,a,t,1] = temp$estimate
      NP_ta_varying[1,a,t,2] = temp$se
      NP_IF_ta_varying[1,a,t,,] = temp$IF
    }
  }
  
  # tg varying
  NP_tg_varying = array(NA, dim = c(2,num_G,J,2))
  NP_IF_tg_varying = array(NA, dim = c(2,num_G,J,n,6))
  dimnames(NP_tg_varying)[[1]] = dimnames(NP_IF_tg_varying)[[1]] = c("et(Xg)","μt(Xg)")
  dimnames(NP_tg_varying)[[2]] = dimnames(NP_IF_tg_varying)[[2]] = G_label
  dimnames(NP_tg_varying)[[3]] = dimnames(NP_IF_tg_varying)[[3]] = T_label
  for (g in 1:num_G) {
    for (t in 1:J) {
      # et(Xg )
      temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                      NP_ti_varying[1,,t],
                      -ones,
                      sampling_weights,cl)
      NP_tg_varying[1,g,t,1] = temp$estimate
      NP_tg_varying[1,g,t,2] = temp$se
      NP_IF_tg_varying[1,g,t,,] = temp$IF
      # μt(Xg )
      temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                      NP_ti_varying[3,,t] #mhat
                      + NP_ti_varying[1,,t] * (Y - NP_ti_varying[3,,t]) / NP_ti_varying[2,,t],
                      -ones,
                      sampling_weights,cl)
      NP_tg_varying[2,g,t,1] = temp$estimate
      NP_tg_varying[2,g,t,2] = temp$se
      NP_IF_tg_varying[2,g,t,,] = temp$IF
    }
  }
  
  # ag varying
  NP_ag_varying = array(NA, dim = c(3,num_A,num_G,2))
  NP_IF_ag_varying = array(NA, dim = c(3,num_A,num_G,n,6))
  dimnames(NP_ag_varying)[[1]] = dimnames(NP_IF_ag_varying)[[1]] = c("eag","ea(Xg)","ma(Xg)")
  dimnames(NP_ag_varying)[[2]] = dimnames(NP_IF_ag_varying)[[2]] = A_label
  dimnames(NP_ag_varying)[[3]] = dimnames(NP_IF_ag_varying)[[3]] = G_label
  for (a in 1:num_A) {
    for (g in 1:num_G) {
      # e_dg
      temp = IF_maker(ones,
                      NP_agi_varying[1,a,,g],
                      -ones,
                      sampling_weights,cl)
      NP_ag_varying[1,a,g,1] = temp$estimate
      NP_ag_varying[1,a,g,2] = temp$se
      NP_IF_ag_varying[1,a,g,,] = temp$IF
      # e_d(g)
      temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                      NP_ai_varying[1,,a],
                      -ones,
                      sampling_weights,cl)
      NP_ag_varying[2,a,g,1] = temp$estimate
      NP_ag_varying[2,a,g,2] = temp$se
      NP_IF_ag_varying[2,a,g,,] = temp$IF
      # m(Xg )
      temp = IF_maker(NP_agi_varying[1,a,,g] / NP_ag_varying[1,a,g,1],
                      Y,
                      -ones,
                      sampling_weights,cl)
      NP_ag_varying[3,a,g,1] = temp$estimate
      NP_ag_varying[3,a,g,2] = temp$se
      NP_IF_ag_varying[3,a,g,,] = temp$IF
    }
  }
  
  # tag varying
  NP_tag_varying = array(NA, dim = c(2,num_A,num_G,J,2))
  NP_IF_tag_varying = array(NA, dim = c(2,num_A,num_G,J,n,6))
  dimnames(NP_tag_varying)[[1]] = dimnames(NP_IF_tag_varying)[[1]] = c("eta(Xg)","E[eta(X)|X ∈ Xg]")
  dimnames(NP_tag_varying)[[2]] = dimnames(NP_IF_tag_varying)[[2]] = A_label
  dimnames(NP_tag_varying)[[3]] = dimnames(NP_IF_tag_varying)[[3]] = G_label
  dimnames(NP_tag_varying)[[4]] = dimnames(NP_IF_tag_varying)[[4]] = T_label
  for (a in 1:num_A) {
    for (g in 1:num_G) {
      for (t in 1:J) {
        # e_td(g)
        temp = IF_maker(NP_agi_varying[1,a,,g] / NP_ag_varying[1,a,g,1],
                        NP_ti_varying[1,,t],
                        -ones,
                        sampling_weights,cl)
        NP_tag_varying[1,a,g,t,1] = temp$estimate
        NP_tag_varying[1,a,g,t,2] = temp$se
        NP_IF_tag_varying[1,a,g,t,,] = temp$IF
        
        # E[eta(X)|X ∈ Xg] 
        temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                        NP_ai_varying[1,,a] * (NP_ti_varying[1,,t] - NP_tai_varying[1,a,,t]) 
                        / NP_ai_varying[2,,a] + NP_tai_varying[1,a,,t],
                        -ones,
                        sampling_weights,cl)
        NP_tag_varying[2,a,g,t,1] = temp$estimate
        NP_tag_varying[2,a,g,t,2] = temp$se
        NP_IF_tag_varying[2,a,g,t,,] = temp$IF
      } # end J
    } # end num_G
  } # end num_A
  
  #########################################################
  ### t-specific components of decomposition parameters
  # tag-varying
  C_tag_varying = array(NA, dim = c(6,num_A,num_G,J,2))
  C_IF_tag_varying = array(NA, dim = c(6,num_A,num_G,J,n,6))
  dimnames(C_tag_varying)[[1]] = dimnames(C_IF_tag_varying)[[1]] = c("eta(Xg)μt(Xg)","E[eta(X)|X ∈ Xg]μt(Xg)","E[eta(X)μt(X)|X ∈ Xg]",
                                                                     "etaμt(Xg)","eta(Xg)μt","Cov(eta(X), μt(X)|X ∈ Xg)")
  dimnames(C_tag_varying)[[2]] = dimnames(C_IF_tag_varying)[[2]] = A_label
  dimnames(C_tag_varying)[[3]] = dimnames(C_IF_tag_varying)[[3]] = G_label
  dimnames(C_tag_varying)[[4]] = dimnames(C_IF_tag_varying)[[4]] = T_label
  for (a in 1:num_A) {
    for (g in 1:num_G) {
      for (t in 1:J) {
        # eta(Xg)μt(Xg)
        temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                        (NP_ai_varying[1,,a] * NP_ti_varying[1,,t]) / NP_ag_varying[2,a,g,1] * NP_tg_varying[2,g,t,1]
                        + NP_tag_varying[1,a,g,t,1] * NP_IF_t_varying[2,t,,5],
                        - (ones + NP_ai_varying[1,,a] / NP_ag_varying[2,a,g,1]),
                        sampling_weights,cl)
        C_tag_varying[1,a,g,t,1] = temp$estimate
        C_tag_varying[1,a,g,t,2] = temp$se
        C_IF_tag_varying[1,a,g,t,,] = temp$IF
        # E[eta(X)|X ∈ Xg]μt(Xg )
        temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                        NP_IF_tag_varying[2,a,g,t,,5] * NP_tg_varying[2,g,t,1]
                        + NP_tag_varying[2,a,g,t,1] * NP_IF_t_varying[2,t,,5],
                        - 2 * ones,
                        sampling_weights,cl)
        C_tag_varying[2,a,g,t,1] = temp$estimate
        C_tag_varying[2,a,g,t,2] = temp$se
        C_IF_tag_varying[2,a,g,t,,] = temp$IF
        # E[eta(X)μt(X)|X ∈ Xg ]
        temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                        NP_ai_varying[1,,a] * (NP_ti_varying[1,,t] - NP_tai_varying[1,a,,t]) 
                        / NP_ai_varying[2,,a] * NP_ti_varying[3,,t]
                        + NP_tai_varying[1,a,,t] * NP_IF_t_varying[2,t,,5],
                        - ones,
                        sampling_weights,cl)
        C_tag_varying[3,a,g,t,1] = temp$estimate
        C_tag_varying[3,a,g,t,2] = temp$se
        C_IF_tag_varying[3,a,g,t,,] = temp$IF
        # etaμt(Xg) 
        temp = IF_maker(ones,
                        NP_ai_varying[1,,a] * NP_ti_varying[1,,t] * NP_tg_varying[2,g,t,1] / NP_a_varying[1,a,1] + 
                          NP_ta_varying[1,a,t,1] * NP_gi_varying[1,,g] / NP_g_varying[1,g,1] * NP_IF_t_varying[2,t,,5],
                        - (NP_ai_varying[1,,a] / NP_a_varying[1,a,1] + NP_gi_varying[1,,g] / NP_g_varying[1,g,1]),
                        sampling_weights,cl)
        C_tag_varying[4,a,g,t,1] = temp$estimate
        C_tag_varying[4,a,g,t,2] = temp$se
        C_IF_tag_varying[4,a,g,t,,] = temp$IF
        # eta(Xg )μt
        temp = IF_maker(ones,
                        NP_agi_varying[1,a,,g] / NP_ag_varying[1,a,g,1]  * NP_ti_varying[1,,t] * NP_t_varying[2,t,1]
                        + NP_tag_varying[1,a,g,t,1] * NP_IF_t_varying[2,t,,5],
                        - (NP_agi_varying[1,a,,g] / NP_ag_varying[1,a,g,1] + ones),
                        sampling_weights,cl)
        C_tag_varying[5,a,g,t,1] = temp$estimate
        C_tag_varying[5,a,g,t,2] = temp$se
        C_IF_tag_varying[5,a,g,t,,] = temp$IF
        # Cov(eta(X), μt(X)|X ∈ Xg )
        C_tag_varying[6,a,g,t,1] = C_tag_varying[3,a,g,t,1] - C_tag_varying[2,a,g,t,1]
        C_IF_tag_varying[6,a,g,t,,1] = C_IF_tag_varying[3,a,g,t,,1] - C_IF_tag_varying[2,a,g,t,,1]
        C_tag_varying[6,a,g,t,2] = IF2SE(C_IF_tag_varying[6,a,g,t,,1],cl)
      }
    }
  }
  
  # tg-varying
  C_tg_varying = array(NA, dim = c(3,num_G,J,2))
  C_IF_tg_varying = array(NA, dim = c(3,num_G,J,n,6))
  dimnames(C_tg_varying)[[1]] = dimnames(C_IF_tg_varying)[[1]] = c("t(Xg)μt(Xg)","E[et(X)μt(X)|X ∈ Xg]","Cov(et(X), μt(X)|X ∈ Xg)")
  dimnames(C_tg_varying)[[2]] = dimnames(C_IF_tg_varying)[[2]] = G_label
  dimnames(C_tg_varying)[[3]] = dimnames(C_IF_tg_varying)[[3]] = T_label
  for (g in 1:num_G) {
    for (t in 1:J) {
      # et(Xg )μt(Xg )
      temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                      NP_ti_varying[1,,t] * NP_tg_varying[2,g,t,1] + NP_tg_varying[1,g,t,1] * NP_IF_t_varying[2,t,,5],
                      - 2 * ones,
                      sampling_weights,cl)
      C_tg_varying[1,g,t,1] = temp$estimate
      C_tg_varying[1,g,t,2] = temp$se
      C_IF_tg_varying[1,g,t,,] = temp$IF
      # E[e(X)μ(X)|X ∈ Xg ]
      temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                      NP_ti_varying[1,,t] * Y,
                      - ones,
                      sampling_weights,cl)
      C_tg_varying[2,g,t,1] = temp$estimate
      C_tg_varying[2,g,t,2] = temp$se
      C_IF_tg_varying[2,g,t,,] = temp$IF
      # Cov(et(X), μt(X)|X ∈ Xg )
      C_tg_varying[3,g,t,1] = C_tg_varying[2,g,t,1] - C_tg_varying[1,g,t,1]
      C_IF_tg_varying[3,g,t,,1] = C_IF_tg_varying[2,g,t,,1] - C_IF_tg_varying[1,g,t,,1]
      C_tg_varying[3,g,t,2] = IF2SE(C_IF_tg_varying[3,g,t,,1],cl)
    }
  }
  
  # ta-varying
  C_ta_varying = array(NA, dim = c(1,num_A,J,2))
  C_IF_ta_varying = array(NA, dim = c(1,num_A,J,n,6))
  dimnames(C_ta_varying)[[1]] = dimnames(C_IF_ta_varying)[[1]] = "etaμt"
  dimnames(C_ta_varying)[[2]] = dimnames(C_IF_ta_varying)[[2]] = A_label
  dimnames(C_ta_varying)[[3]] = dimnames(C_IF_ta_varying)[[3]] = T_label
  for (a in 1:num_A) {
    for (t in 1:J) {
      # etaμt
      temp = IF_maker(ones,
                      NP_ai_varying[1,,a] / NP_a_varying[1,a,1] * NP_ti_varying[1,,t] * NP_t_varying[2,t,1] + 
                        NP_ta_varying[1,a,t,1] * NP_IF_t_varying[2,t,,5] ,
                      - (NP_ai_varying[1,,a] / NP_a_varying[1,a,1] + 1),
                      sampling_weights,cl)
      C_ta_varying[1,a,t,1] = temp$estimate
      C_ta_varying[1,a,t,2] = temp$se
      C_IF_ta_varying[1,a,t,,] = temp$IF
    }
  }
  
  
  #########################################################
  ### Decomposition target parameters - levels
  TP_level= array(NA, dim = c(14,num_G,num_A,2))
  TP_level_IF = array(NA, dim = c(14,num_G,num_A,n,6))
  dimnames(TP_level)[[1]] = dimnames(TP_level_IF)[[1]] = c("GM","AGM","d0","s1","s2","s3","d1","d2","d3",
                                                           "Cov(etX,mut|Xg)","d4","SRCT2","d4'","d5")
  dimnames(TP_level)[[2]] = dimnames(TP_level_IF)[[2]] = G_label
  dimnames(TP_level)[[3]] = dimnames(TP_level_IF)[[3]] = A_label
  dimnames(TP_level)[[4]] = c("PE","SE")
  for (g in 1:num_G) {
    for (a in 1:num_A) {
      ### Starting points to be decomposed
      # ma(Xg)
      temp = IF_maker(NP_agi_varying[1,a,,g] / NP_ag_varying[1,a,g,1],
                      Y,
                      -ones,
                      sampling_weights,cl)
      TP_level[1,g,a,1] = temp$estimate
      TP_level[1,g,a,2] = temp$se
      TP_level_IF[1,g,a,,] = temp$IF
      # μa(Xg)
      temp = IF_maker(NP_gi_varying[1,,g] / NP_g_varying[1,g,1],
                      NP_ai_varying[1,,a] * (Y - NP_ai_varying[3,,a]) / NP_ai_varying[2,,a] + NP_ai_varying[3,,a],
                      -ones,
                      sampling_weights,cl)
      TP_level[2,g,a,1] = temp$estimate
      TP_level[2,g,a,2] = temp$se
      TP_level_IF[2,g,a,,] = temp$IF
      
      ### Decomposition components
      # d0: eta*μt
      temp = IF_maker(ones, 
                      colSums(C_IF_ta_varying[1,a,,,3]),
                      colMeans(C_IF_ta_varying[1,a,,,2]),
                      sampling_weights,cl)
      TP_level[3,g,a,1] = temp$estimate
      TP_level[3,g,a,2] = temp$se
      TP_level_IF[3,g,a,,] = temp$IF
      
      # eta*μt(Xg)
      temp = IF_maker(ones, 
                      colSums(C_IF_tag_varying[4,a,g,,,3]),
                      colMeans(C_IF_tag_varying[4,a,g,,,2]),
                      sampling_weights,cl)
      TP_level[4,g,a,1] = temp$estimate
      TP_level[4,g,a,2] = temp$se
      TP_level_IF[4,g,a,,] = temp$IF
      
      # eta(Xg)*μt
      temp = IF_maker(ones, 
                      colSums(C_IF_tag_varying[5,a,g,,,3]),
                      colMeans(C_IF_tag_varying[5,a,g,,,2]),
                      sampling_weights,cl)
      TP_level[5,g,a,1] = temp$estimate
      TP_level[5,g,a,2] = temp$se
      TP_level_IF[5,g,a,,] = temp$IF
      
      # eta(Xg)*μt(Xg)
      temp = IF_maker(ones, 
                      colSums(C_IF_tag_varying[1,a,g,,,3]),
                      colMeans(C_IF_tag_varying[1,a,g,,,2]),
                      sampling_weights,cl)
      TP_level[6,g,a,1] = temp$estimate
      TP_level[6,g,a,2] = temp$se
      TP_level_IF[6,g,a,,] = temp$IF
      
      ### Decomposition parameters
      # d1 
      TP_level[7,g,a,1] = TP_level[4,g,a,1] - TP_level[3,g,a,1]
      TP_level_IF[7,g,a,,1] = TP_level_IF[4,g,a,,1] - TP_level_IF[3,g,a,,1]
      TP_level[7,g,a,2] = IF2SE(TP_level_IF[7,g,a,,1],cl)
      
      # d2
      TP_level[8,g,a,1] = TP_level[5,g,a,1] - TP_level[3,g,a,1]
      TP_level_IF[8,g,a,,1] = TP_level_IF[5,g,a,,1] - TP_level_IF[3,g,a,,1]
      TP_level[8,g,a,2] = IF2SE(TP_level_IF[8,g,a,,1],cl)
      
      # d3
      TP_level[9,g,a,1] = TP_level[6,g,a,1] - TP_level[4,g,a,1] - TP_level[5,g,a,1] + TP_level[3,g,a,1]
      TP_level_IF[9,g,a,,1] = TP_level_IF[6,g,a,,1] - TP_level_IF[4,g,a,,1] - TP_level_IF[5,g,a,,1] + TP_level_IF[3,g,a,,1]
      TP_level[9,g,a,2] = IF2SE(TP_level_IF[9,g,a,,1],cl)
      
      # Cov(etX,mutX|XG)
      TP_level[10,g,a,1] = sum(C_tg_varying[3,g,treatment_mapping[,a],1])
      # Handle common special case where column dimension gets lost and drop=F would keep all array dims
      temp = C_IF_tg_varying[3,g,treatment_mapping[,a],,1]
      if (is.null(dim(temp))) TP_level_IF[10,g,a,,1] = temp
      else TP_level_IF[10,g,a,,1] = colSums(temp)
      TP_level[10,g,a,2] = IF2SE(TP_level_IF[10,g,a,,1],cl)
      
      # d4: Cov(etX,mutX|XG) / ea(Xg)
      TP_level[11,g,a,1] = TP_level[10,g,a,1] / NP_ag_varying[2,a,g,1]
      TP_level_IF[11,g,a,,1] = 1 / NP_ag_varying[2,a,g,1] * TP_level_IF[10,g,a,,1] - 
        TP_level[10,g,a,1] / NP_ag_varying[2,a,g,1]^2 * NP_IF_ag_varying[2,a,g,,1]
      TP_level[11,g,a,2] = IF2SE(TP_level_IF[11,g,a,,1],cl)
      
      # E[eta(X)|X ∈ Xg]μt(Xg )
      temp = IF_maker(ones, 
                      colSums(C_IF_tag_varying[2,a,g,,,3]),
                      colMeans(C_IF_tag_varying[2,a,g,,,2]),
                      sampling_weights,cl)
      TP_level[12,g,a,1] = temp$estimate
      TP_level[12,g,a,2] = temp$se
      TP_level_IF[12,g,a,,] = temp$IF
    
      # d4': Cov(etaX,mutX|XG)
      TP_level[13,g,a,1] = sum(C_tag_varying[6,a,g,treatment_mapping[,a],1])
      # Handle special case where column dimension gets lost and drop=F would keep all array dims
      temp = C_IF_tag_varying[6,a,g,treatment_mapping[,a],,1]
      if (is.null(dim(temp))) TP_level_IF[13,g,a,,1] = temp
      else TP_level_IF[13,g,a,,1] = colSums(temp)
      TP_level[13,g,a,2] = IF2SE(TP_level_IF[13,g,a,,1],cl)
      
      # d5: E[eta(X)|X ∈ Xg]μt(Xg ) - eta(Xg )μt(Xg )
      TP_level[14,g,a,1] = TP_level[12,g,a,1] - TP_level[6,g,a,1]
      TP_level_IF[14,g,a,,1] = TP_level_IF[12,g,a,,1] - TP_level_IF[6,g,a,,1]
      TP_level[14,g,a,2] = IF2SE(TP_level_IF[14,g,a,,1],cl)
    }
  }
  # Sanity check that decomposition parameters add up to starting points
  # TP_level[3,,,1] + TP_level[7,,,1] + TP_level[8,,,1] + TP_level[9,,,1] + TP_level[11,,,1]
  # TP_level[1,,,1]
  # TP_level[3,,,1] + TP_level[7,,,1] + TP_level[8,,,1] + TP_level[9,,,1] + TP_level[13,,,1] + TP_level[14,,,1]
  # TP_level[2,,,1]
  
  output = list(parameter = TP_level,
                IFs = TP_level_IF,
                mapping = treatment_mapping, 
                label = list(A_label,T_label,G_label),
                cl = cl)
  
  class(output) = "HK2_decomposition"
  output
}



#####################################################
###### Utils

#' Finds mapping between aggregated and effective treatment
#'
#' @param A_mat Logical matrix of aggregated treatment. 
#' For example created by \code{\link{prep_w_mat}}.
#' @param T_mat Logical matrix of effective treatment.
#' For example created by \code{\link{prep_w_mat}}.
#'
#' @return Logical matrix storing the mapping of aggregate and effective treatment.
#'
#' @keywords internal
#' @noRd
#'
treatment_aggregation_mapping = function(A_mat, T_mat) {
  mapping = t(T_mat) %*% A_mat > 0
  if (sum( rowSums(mapping) > 1 ) > 0 ) stop("Provided treatment aggregation not mutually exclusive.")
  return(mapping)
}


#' Statistical inference for parameter with influence functions of the form
#' Psi = θ Psia + Psib = PsiZ (PsiY + θ PsiD).
#'
#' @param PsiZ PsiZ component.
#' @param PsiY PsiY component.
#' @param PsiD PsiD component.
#' @param weights Optional vector of sampling weights.
#' @param cl If not NULL, vector with cluster variables.
#'
#' @return List of three components:
#' - \code{estimate} Point estimate.
#' - \code{se} Standard error.
#' - \code{IF} n x 6 matrix with columns Psi, Psia, Psib, PsiZ, PsiY, PsiD.
#'
#' @keywords internal
#' @noRd
#'
IF_maker = function(PsiZ,PsiY,PsiD,
                    weights=NULL,
                    cl = NULL) {
  N = length(PsiY)
  
  if (is.null(weights)) weights = rep(1,N)
  
  weights = weights / sum(weights) * N
  
  Psia = PsiD * PsiZ * weights
  Psib = PsiY * PsiZ * weights
  estimate = -sum(Psib) / sum(Psia)
  Psi = estimate * Psia + Psib
  se = IF2SE(Psi,cl)
  IF = cbind(Psi,Psia,Psib,PsiZ,PsiY,PsiD)
  colnames(IF) = c("Psi","Psia","Psib","PsiZ","PsiY","PsiD")
  output = list("estimate" = unlist(estimate), "se" = unlist(se), "IF" = IF)
  return(output)
}


#' Calculates either plain or cluster robust standard error from IF.
#'
#' @param IF Influence function vector
#' @param cl If not NULL, vector with cluster variables.
#'
#' @return Standard error
#'
#' @keywords internal
#' @noRd
#'
IF2SE = function(IF, cl = NULL) {
  N = length(IF)
  
  if(is.null(cl)) {
    # Standard standard error.
    se = sqrt(sum(IF^2)) / N
  } else {
    if(length(cl) != N) stop("'cl' must be the same length as the influence function vector 'IF'.")
    # Sum of the influence function values within each cluster.
    cluster_sum = tapply(IF, cl, sum)
    # Cluster-robust standard error.
    se = sqrt(sum(cluster_sum^2)) / N
  }
  return(se)
}



#' Tailored collection of parameters produced by \code{\link{HK2_decomposition}}.
#'
#' @param HK2_decomposition Object of class \code{\link{HK2_decomposition}}.
#' @param parameter Parameter to be collected. Either "adim" or "dim" for the respective decomposition 
#' or a subset of 1:14. Default "adim".
#' @param t_aggregate Integer scalar of vector of length two specifying treatment aggregate to be collected. 
#' 1 is reserved for collection of all treaments. 
#' If single integer, level is extracted. 
#' Vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#' @param x_aggregate Integer scalar of vector of length two specifying covariate aggregate to be collected. 
#' If single integer, level is extracted. 
#' 1 is reserved for the unconditional version. 
#' vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#'
#' @export
#'
HK2_result_collector = function(HK2_decomposition, 
                            parameter = "adim",
                            t_aggregate = c(3,2),
                            x_aggregate = c(3,2)) {
  if (all(parameter %in% 1:14)) params = parameter
  else if (parameter == "dim") params = c(1,3,7,8,9,11)
  else if (parameter == "adim") params = c(2,3,7,8,9,13,14)
  else stop("Please provide valid parameter option.")
  
  res = matrix(NA,nrow = length(params), ncol = 4)
  colnames(res) = c("Estimate","SE","t-value","p-value")
  temp = HK2_parameter_maker(HK2_decomposition$parameter,
                            HK2_decomposition$IFs,
                            parameter = params,
                            t_aggregate = t_aggregate,
                            x_aggregate = x_aggregate,
                            cl = HK2_decomposition$cl)
  res[,1:2] = temp
  rownames(res) = rownames(temp)
  # To immunize against machine noise huge t-values
  res[abs(res) < 1.0e-10] = 0
  
  # t-stat
  res[,3] = res[,1] / res[,2]
  res[,3][is.na(res[,3])] = 0
  # p-value
  res[,4] = 2 * stats::pt(abs(res[,3]),nrow(HK2_decomposition$IFs),lower = FALSE)
  return(res)
}


#' \code{summary} method for class \code{\link{HK2_decomposition}}
#'
#' @param HK2_decomposition Object of class \code{\link{HK2_decomposition}}.
#' @param parameter Parameter to be collected. Either "adim" or "dim" for the respective decomposition 
#' or a subset of 1:14. Default "adim".
#' @param t_aggregate Integer scalar of vector of length two specifying treatment aggregate to be collected. 
#' 1 is reserved for collection of all treaments. 
#' If single integer, level is extracted. 
#' Vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#' @param x_aggregate Integer scalar of vector of length two specifying covariate aggregate to be collected. 
#' If single integer, level is extracted. 
#' 1 is reserved for the unconditional version. 
#' vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#'
#' @return Matrix containing point estimates and standard errors of selected parameters.
#'
#' @export
#'
summary.HK2_decomposition = function(HK2_decomposition, 
                                  parameter = "adim",
                                  t_aggregate = c(2,1),
                                  x_aggregate = c(3,2)) {
  res = HK2_result_collector(HK2_decomposition, 
                         parameter = parameter,
                         t_aggregate = t_aggregate,
                         x_aggregate = x_aggregate)
  
  HK2_message_maker(HK2_decomposition$label,
                   parameter = parameter,
                   t_aggregate = t_aggregate,
                   x_aggregate = x_aggregate)
  
  printCoefmat(res,has.Pvalue = TRUE)
  
  return(invisible(res))
}


#' \code{print} method for class \code{\link{HK2_decomposition}} summarizing the setting.
#'
#' @param HK2_decomposition Object of class \code{\link{HK2_decomposition}}.
#' 
#' @export
#'
print.HK2_decomposition = function(HK2_decomposition) {
  # Create the flexible print statement using sprintf
  message = sprintf("================= Setting Summary =================\n")
  # Print the message
  cat(message)
  for (i in 2:length(HK2_decomposition$label[[1]])) {
    message = sprintf("Effective treatments: %s\n", 
                      paste(HK2_decomposition$label[[2]][HK2_decomposition$mapping[,i]], collapse = ", "))
    cat(message)
    message = sprintf("Aggregated into: %s\n\n", 
                      HK2_decomposition$label[[1]][i])
    cat(message)
  }
  message = sprintf("Heterogeneity variables: %s\n\n", paste(HK2_decomposition$label[[3]], collapse = ", "))
  cat(message)
}


#' \code{plot} method for class \code{\link{HK2_decomposition}}
#'
#' @param HK2_decomposition Object of class \code{\link{HK2_decomposition}}.
#' @param decomposition Either "adim" or "dim" for the respective decomposition. Default "adim".
#' @param t_aggregate Integer scalar of vector of length two specifying treatment aggregate to be collected. 
#' 1 is reserved for collection of all treaments. 
#' If single integer, level is extracted. 
#' Vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#' @param x_aggregate Integer scalar of vector of length two specifying covariate aggregate to be collected. 
#' If single integer, level is extracted. 
#' 1 is reserved for the unconditional version. 
#' vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#' @param levels If TRUE, prints levels next to difference. Only applicable with length(x_aggregate) == 2.
#' @param pe_digits Controls number of digits printed for point estimate.
#' @param distance Controls distance between point estimate and p-value for fine-tuning. Default 0.04.
#' 
#' @importFrom dplyr %>%
#' 
#' @export
#'
plot.HK2_decomposition = function(HK2_decomposition, 
                                  decomposition = "adim",
                                  t_aggregate = c(3,2),
                                  x_aggregate = c(3,2),
                                  levels = F,
                                  pe_digits = 3,
                                  distance = 0.04) {
  
  # Hard-dcode colors obtained from viridis(7,begin = 0.3,alpha=0.9)
  colors = c("#35608DE6", "#287C8EE6", "#1F988BE6", "#2FB47CE6",
             "#66CB5DE6", "#B1DD2FE6", "#FDE725E6")
  
  if (decomposition == "dim") {selector = c(3,7,8,9,11,1); colors = colors[-6]}
  else if (decomposition == "adim") selector = c(3,7,8,9,13,14,2)
  else stop("Please provide valid decomposition option.")
  
  temp = HK2_result_collector(HK2_decomposition,parameter = selector, 
                              t_aggregate = t_aggregate,x_aggregate = x_aggregate)
  # Extract coefficients
  coef = temp[,1]
  pvalue = temp[,4]
  format(pvalue,digits=1,nsmall=3,scientific=F)
  # Format results
  coef = round(coef,pe_digits)
  pvalue = paste0("(",format(round(pvalue,3),digits=1,nsmall=3,scientific=F),")")
  # Manual fix
  pvalue[pvalue == "(0.000)"] = "(<0.001)"
  
  # Start and end points for waterfall graph
  start = c(0, cumsum(coef[-length(coef)]))
  end = c(cumsum(coef[-length(coef)]),0)
  
  estimand_labels = names(coef)
  
  if (isFALSE(levels)) {
    range_values = max(start) - min(start)
    
    # Plot Waterfall
    data.frame(id=1:length(end),
               Estimand = factor(estimand_labels,levels = estimand_labels),
               start = start,
               end = end,
               coef = sprintf("%.1f", coef),
               pvalue = pvalue) %>%
      ggplot(aes(id, fill = Estimand)) +
      geom_rect(aes(xmin = id - 0.5, xmax = id + 0.5, ymin = end,ymax = start)) + 
      scale_fill_manual(values = colors) + ylab("Decomposition") + 
      scale_x_continuous(breaks=1:length(estimand_labels),labels=estimand_labels) + 
      theme_bw() + geom_hline(yintercept = 0) +
      xlab("Estimand") + theme(legend.position="none") +
      geom_text(aes(x = id, y = (start+end+range_values*distance)/2, label = coef),size = 2.75) +
      geom_text(aes(x = id, y = (start+end-range_values*distance)/2, label = pvalue),size = 2.75)
  }
  else if (isTRUE(levels) & length(x_aggregate) == 2) {
    # First level
    temp = HK2_result_collector(HK2_decomposition,parameter = selector, 
                                t_aggregate = t_aggregate,x_aggregate = x_aggregate[1])
    # Extract coefficients
    coef1 = temp[,1]
    pvalue1 = temp[,4]
    format(pvalue1,digits=1,nsmall=3,scientific=F)
    # Format results
    coef1 = round(coef1,pe_digits)
    pvalue1 = paste0("(",format(round(pvalue1,3),digits=1,nsmall=3,scientific=F),")")
    # Manual fix
    pvalue1[pvalue1 == "(0.000)"] = "(<0.001)"
    
    # Start and end points for waterfall graph
    start1 = c(0, cumsum(coef1[-length(coef1)]))
    end1 = c(cumsum(coef1[-length(coef1)]),0)
    
    # Second level
    temp = HK2_result_collector(HK2_decomposition,parameter = selector, 
                                t_aggregate = t_aggregate,x_aggregate = x_aggregate[2])
    # Extract coefficients
    coef2 = temp[,1]
    pvalue2 = temp[,4]
    format(pvalue2,digits=1,nsmall=3,scientific=F)
    # Format results
    coef2 = round(coef2,pe_digits)
    pvalue2 = paste0("(",format(round(pvalue2,3),digits=1,nsmall=3,scientific=F),")")
    # Manual fix
    pvalue2[pvalue2 == "(0.000)"] = "(<0.001)"
    
    # Start and end points for waterfall graph
    start2 = c(0, cumsum(coef2[-length(coef2)]))
    end2 = c(cumsum(coef2[-length(coef2)]),0)
    
    label1 = HK2_decomposition$label[[3]][x_aggregate[1]]
    label2 = HK2_decomposition$label[[3]][x_aggregate[2]]
    label = paste("Contrast",label1,"-",label2)
    
    range_values = max(c(start1,start2,start)) - min(c(start1,start2,start))
    
    # Plot Waterfall
    # Customized x-axis require a bit of tweaking
    labels_group1 = names(coef1)
    labels_group2 = names(coef2)
    labels_group3 = names(coef)
    
    # Create a new variable 'xlab' that depends on Group.
    data = data.frame(
      id = rep(1:length(end), 3),
      Estimand = factor(rep(estimand_labels, 3), levels = estimand_labels),
      Group = factor(c(rep(label1, length(end)), rep(label2, length(end)), rep(label, length(end))),
                     levels = c(label1, label2, label)),
      start = c(start1, start2, start),
      end = c(end1, end2, end),
      coef = sprintf("%.1f", c(coef1, coef2, coef)),
      pvalue = c(pvalue1, pvalue2, pvalue)
    )
    unique_labels = c(labels_group1[1:length(coef1)-1],labels_group3)
    # Create a facet-specific x-axis label variable.
    data$xlab = factor(with(data, ifelse(Group == label1, 
                                         labels_group1[id],
                                         ifelse(Group == label2,
                                                labels_group2[id],
                                                labels_group3[id]))),
                       levels = unique_labels)
    
    # Now plot using the discrete x variable and free x scales.
    ggplot(data, aes(x = xlab, fill = Estimand, group = Group)) +
      geom_rect(aes(xmin = as.numeric(id) - 0.5, xmax = as.numeric(id) + 0.5, 
                    ymin = end, ymax = start)) +
      scale_fill_manual(values = colors) +
      ylab("Decomposition") +
      theme_bw() +
      geom_hline(yintercept = 0) +
      xlab("Estimand") +
      theme(legend.position = "none") +
      geom_text(aes(x = xlab, y = (start + end + range_values * distance) / 2, label = coef), size = 2.75) +
      geom_text(aes(x = xlab, y = (start + end - range_values * distance) / 2, label = pvalue), size = 2.75) +
      facet_wrap(~Group, scales = "free_x")
    
  }
  else stop("levels = TRUE only applicable with length(x_aggregate) == 2")
}



#' Internal function to calculate point estimates and IF based
#' standard errors from the ag-specific levels.
#'
#' @param level Level point estimates.
#' @param parameter Parameter to be collected. Either "adim" or "dim" for the respective decomposition 
#' or a subset of 1:14. Default "adim".
#' @param t_aggregate Integer scalar of vector of length two specifying treatment aggregate to be collected. 
#' 1 is reserved for collection of all treaments. 
#' If single integer, level is extracted. 
#' Vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#' @param x_aggregate Integer scalar of vector of length two specifying covariate aggregate to be collected. 
#' If single integer, level is extracted. 
#' 1 is reserved for the unconditional version. 
#' vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#' @param cl If not NULL, vector with cluster variables.
#'
#' @return Matrix containing point estimates and standard errors of selected parameters.
#'
#' @keywords internal
#' @noRd
#'
HK2_parameter_maker = function(level,
                               level_IF,
                               parameter = 1:14,
                               t_aggregate,
                               x_aggregate,
                               cl = NULL) {
  results = matrix(NA,nrow = 14, ncol = 2)
  # IFs = matrix(NA,nrow = nrow(level_IF[1,1,1,,]), ncol = 14) # for sanity check below
  colnames(results) = c("PE","SE")
  for (p in 1:14) {
    if (length(t_aggregate) == 1 & length(x_aggregate) == 1) {
      results[p,1] = level[p,x_aggregate,t_aggregate,1]
      IF = level_IF[p,x_aggregate,t_aggregate,,]
      rownames(results) = c("GM","AGM","d0","s1","s2","s3","d1","d2","d3",
                            "Cov(etX,mut|Xg)","d4","SRCT2","d4'","d5")
    }
    if (length(t_aggregate) == 2 & length(x_aggregate) == 1) {
      results[p,1] = level[p,x_aggregate,t_aggregate[1],1] - level[p,x_aggregate,t_aggregate[2],1]
      IF = level_IF[p,x_aggregate,t_aggregate[1],,] - level_IF[p,x_aggregate,t_aggregate[2],,]
      rownames(results) = c("DiM","ADiM","δ0","s1","s2","s3","δ1","δ2","δ3",
                            "Cov(etX,mut|Xg)","δ4","SRCT2","δ4'","δ5")
    }
    if (length(t_aggregate) == 1 & length(x_aggregate) == 2) {
      results[p,1] = level[p,x_aggregate[1],t_aggregate,1] - level[p,x_aggregate[2],t_aggregate,1]
      IF = level_IF[p,x_aggregate[1],t_aggregate,,] - level_IF[p,x_aggregate[2],t_aggregate,,]
      rownames(results) = c("GM","AGM","d0","l1","l2","l3","d1","d2","d3",
                            "Cov(etX,mut|Xg)","d4","SRCT2","d4'","d5")
    }
    if (length(t_aggregate) == 2 & length(x_aggregate) == 2) {
      results[p,1] = level[p,x_aggregate[1],t_aggregate[1],1] - level[p,x_aggregate[1],t_aggregate[2],1] - 
        (level[p,x_aggregate[2],t_aggregate[1],1] - level[p,x_aggregate[2],t_aggregate[2],1])
      IF = level_IF[p,x_aggregate[1],t_aggregate[1],,] - level_IF[p,x_aggregate[1],t_aggregate[2],,] - 
        (level_IF[p,x_aggregate[2],t_aggregate[1],,] - level_IF[p,x_aggregate[2],t_aggregate[2],,])
      rownames(results) = c("DiM","ADiM","Δ0","l1","l2","l3","Δ1","Δ2","Δ3",
                            "Cov(etX,mut|Xg)","Δ4","SRCT2","Δ4'","Δ5")
    }
    results[p,2] = IF2SE(IF[,1],cl)
    # IFs[,p] = IF[,1]
  }
  # Sanity check that decomposition parameters add up to starting points including same SEs
  # print(sum(results[c(3,7,8,9,11),1]))
  # print(sum(results[c(3,7,8,9,13,14),1]))
  # print( sqrt( var(rowSums(IFs[,c(3,7,8,9,11)])) / nrow(IF) ) )
  # print( sqrt( var(rowSums(IFs[,c(3,7,8,9,13,14)])) / nrow(IF) ) )
  
  return(results[parameter,])
}


#' Internal function to print the setting for which results are display.
#'
#' @param label A_label, T_label, and G_label.
#' @param parameter Parameters to be collected.
#' @param t_aggregate Integer scalar of vector of length two specifying treatment aggregate to be collected. 
#' 1 is reserved for collection of all treaments. 
#' If single integer, level is extracted. 
#' Vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#' @param x_aggregate Integer scalar of vector of length two specifying covariate aggregate to be collected. 
#' If single integer, level is extracted. 
#' 1 is reserved for the unconditional version. 
#' vector of length two like c(3,2) provides aggregate 2 minus aggregate 1.
#'
#' @return None. Just prints a message.
#'
#' @keywords internal
#' @noRd
#'
HK2_message_maker = function(label,
                             parameter = 1:14,
                             t_aggregate,
                             x_aggregate) {
  message = sprintf("================= HK decomposition results for =================\nParameters %s\n\n", 
                    paste(parameter, collapse = ", "))
  cat(message)
  if (length(t_aggregate) == 1 & length(x_aggregate) == 1) {
    message = sprintf("Treatment aggregation level: %s\n\n", 
                      label[[1]][t_aggregate])
    cat(message)
    message = sprintf("Heterogeneity level: %s\n\n", 
                      label[[3]][x_aggregate])
    cat(message)
  }
  if (length(t_aggregate) == 2 & length(x_aggregate) == 1) {
    message = sprintf("Treatment aggregation contrast: %s - %s\n\n", 
                      label[[1]][t_aggregate[1]],label[[1]][t_aggregate[2]])
    cat(message)
    message = sprintf("Heterogeneity level: %s\n\n", 
                      label[[3]][x_aggregate])
    cat(message)
  }
  if (length(t_aggregate) == 1 & length(x_aggregate) == 2) {
    message = sprintf("Treatment aggregation level: %s\n\n", 
                      label[[1]][t_aggregate])
    cat(message)
    message = sprintf("Heterogeneity contrast: %s - %s \n\n", 
                      label[[3]][x_aggregate[1]],label[[3]][x_aggregate[2]])
    cat(message)
  }
  if (length(t_aggregate) == 2 & length(x_aggregate) == 2) {
    message = sprintf("Treatment aggregation contrast: %s - %s\n\n", 
                      label[[1]][t_aggregate[1]],label[[1]][t_aggregate[2]])
    cat(message)
    message = sprintf("Heterogeneity contrast: %s - %s \n\n", 
                      label[[3]][x_aggregate[1]],label[[3]][x_aggregate[2]])
    cat(message)
  }
}



