#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat updateEta_c(arma::mat Lambda, arma::vec ps, int k, arma::mat Y, int n) {
  mat Lmsg = Lambda.each_col() % ps;
  mat Veta1 = eye<mat>(k, k) + trans(Lambda) * Lmsg;

  // Cholesky decomposition
  mat S = chol(Veta1);

  // Matrix inversion using the Cholesky decomposition
  mat Veta = inv_sympd(Veta1);

  mat Meta = Y * Lmsg * Veta;

  // Generate random matrix
  mat randMat = randn(n, k);

  // Update eta
  mat eta = Meta + randMat * S;

  return eta;
}

// [[Rcpp::export]]
arma::mat updateLambda_c(arma::mat eta, arma::mat Plam, arma::vec ps, arma::mat Y, int k, int p) {
  arma::mat eta2 = eta.t() * eta;
  arma::mat Lambda(p, k, arma::fill::zeros);

  for (int j = 0; j < p; ++j) {
    arma::mat Qlam = diagmat(Plam.row(j)) + ps(j) * eta2;
    arma::mat blam = ps(j) * eta.t() * Y.col(j);
    arma::mat Llam = chol(Qlam);
    arma::mat zlam = arma::randn<arma::mat>(k, 1);
    arma::mat vlam = inv_sympd(Qlam);
    arma::mat mlam = vlam * blam;
    arma::mat ylam = vlam * zlam;
    Lambda.row(j) = (ylam + mlam).t();
  }

  return Lambda;
}

// [[Rcpp::export]]
arma::mat updatePsi_c(double df, arma::mat Lambda, arma::vec tauh, int p, int k) {
  mat Lambda2 = Lambda % Lambda;
  mat df2 = df / 2 + Lambda2.each_row() % tauh.t() / 2;
  mat psijh(p, k, fill::zeros);

  for (int j = 0; j < p; ++j) {
    for (int h = 0; h < k; ++h) {
      psijh(j, h) = R::rgamma(df / 2 + 0.5, 1.0 / df2(j, h));
    }
  }

  return psijh;
}

// [[Rcpp::export]]
List updateDeltaTauh_c(arma::mat Lambda, arma::mat psijh, double ad1, int p, int k,
                       double bd1, arma::vec delta, arma::vec tauh, double ad2, double bd2) {

  mat temp = Lambda % Lambda % psijh;
  vec tmp = sum(temp, 0).t();  // Sum along the rows and transpose
  vec newDelta = delta;
  vec newTauh = tauh;

  double ad, bd;

  ad = ad1 + 0.5 * p * k;
  bd = bd1 + 0.5 * sum(tauh % tmp) / delta(0);
  newDelta(0) = R::rgamma(ad, 1 / bd);
  newTauh = cumprod(newDelta);

  for (int h = 1; h < k; ++h) {
    ad = ad2 + 0.5 * p * (k - h);
    bd = bd2 + 0.5 / delta(h) * sum(newTauh.subvec(h, k - 1) % tmp.subvec(h, k - 1));
    newDelta(h) = R::rgamma(ad, 1 / bd);
    newTauh = cumprod(newDelta);
  }

  return List::create(Named("delta") = newDelta, Named("tauh") = newTauh);
}

// [[Rcpp::export]]
arma::vec updateSigma_c(arma::vec tmp2, int p, double as, int n, double bs) {

  // Initialize ps vector
  arma::vec ps(p);

  // Update each element of ps
  for (int j = 0; j < p; ++j) {
    ps(j) = R::rgamma(as + 0.5 * n, 1 / (bs + 0.5 * tmp2(j)));
  }

  return ps;
}


// adaptation in  slows down
