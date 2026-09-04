#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// The spacetime operator's K is assembled in C++ (Spacetime::build_KZ) and in R
// (the update_K closure in complex-models.R). Both must agree: R builds the
// initial K and drives simulate(), C++ drives every fitting iteration. This
// binding exposes the C++ assembly so a test can compare them directly, rather
// than inferring agreement from whether a fit happens to converge.
#include <latents/tensorprod.cpp>

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::SparseMatrix<double> cpp_spacetime_K(const List &operator_list,
                                            const Eigen::VectorXd &theta_K) {
  Spacetime op(operator_list);
  op.build_KZ(theta_K);
  return op.getK();
}
