#include "hiperglm_types.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List eigenLeastSquaresWeighted(const Eigen::Map<MatrixXd> &X,
                               const Eigen::Map<VectorXd> &y,
                               const Eigen::Map<VectorXd> &w) {
  VectorXd sqrt_w = w.array().sqrt();
  MatrixXd tilde_X = X;
  for (int i = 0; i < X.rows(); i++) {
    tilde_X.row(i) *= sqrt_w(i);
  }
  VectorXd tilde_y = y.array() * sqrt_w.array();
  Eigen::HouseholderQR<MatrixXd> qr(tilde_X);
  VectorXd beta = qr.solve(tilde_y);
  VectorXd fitted = X * beta;
  VectorXd resid = y - fitted;
  return List::create(
    Named("coef") = beta,
    Named("fitted") = fitted,
    Named("resid") = resid
  );
}
