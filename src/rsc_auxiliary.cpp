#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec rsc_d2subspace(arma::mat X, arma::mat projection, arma::vec center){
  // parameters
  int N = X.n_rows;
  int p = X.n_cols;
  
  // prepare
  // arma::mat proj = arma::eye<arma::mat>(p,p) - projection;
  
  arma::mat Ip(p,p,fill::eye);
  arma::mat pr = Ip - projection;
  arma::vec xdiff(p,fill::zeros);
  arma::vec dvec(N,fill::zeros);
  // compute
  for (unsigned int n=0;n<N;n++){
    xdiff = arma::trans(X.row(n)) - center;
    dvec(n) = arma::norm(pr*xdiff, 2);
  }
  
  // return
  return(dvec);
}