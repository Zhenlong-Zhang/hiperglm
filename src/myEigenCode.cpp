#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd eigenTimesTwo(const Eigen::VectorXd &x) {
  return x * 2;
}
