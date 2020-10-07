#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec norm_drl_rcpp(arma::sp_mat alpha, arma::mat y_mat, arma::colvec y, arma::mat w_mat, arma::mat e_mat) {

  int n = alpha.n_rows;
  arma::mat alphad(alpha);
  arma::mat correction = alphad * (w_mat / e_mat);
  // cout << "Maximum correction " << max(correction) << endl;
  arma::mat res = -y_mat;
  res.each_col() += y;
  arma::colvec diff_maker = { -1, 1 };
  arma::mat temp(n,2);
  arma::vec cates_norm(n);

  for(int i = 0; i < n; i++) {
    // cout << "Row " << i << endl;
    arma::mat C = arma::diagmat(correction.row(i));
    temp = y_mat + res % (w_mat / (e_mat * C));
    cates_norm(i) = arma::as_scalar(alphad.row(i) * temp * diff_maker);
    Rcpp::checkUserInterrupt();
  }

  return cates_norm;
}
